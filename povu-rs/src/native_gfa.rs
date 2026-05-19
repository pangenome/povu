use crate::{
    AlleleConstruction, AlternateAllele, Contig, Error, GenotypeAllele, GenotypeColumn, Result,
    VariantCall, VariantSource, VariantType, VcfDocument,
};
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    fn from_gfa(value: &str) -> Result<Self> {
        match value {
            "+" => Ok(Self::Forward),
            "-" => Ok(Self::Reverse),
            other => Err(Error::invalid_gfa(format!(
                "unsupported GFA orientation '{other}'"
            ))),
        }
    }

    fn from_step_suffix(value: char) -> Result<Self> {
        match value {
            '+' => Ok(Self::Forward),
            '-' => Ok(Self::Reverse),
            other => Err(Error::invalid_gfa(format!(
                "unsupported path step orientation '{other}'"
            ))),
        }
    }

    fn flip(self) -> Self {
        match self {
            Self::Forward => Self::Reverse,
            Self::Reverse => Self::Forward,
        }
    }

    fn token_prefix(self) -> char {
        match self {
            Self::Forward => '>',
            Self::Reverse => '<',
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Step {
    pub segment: String,
    pub strand: Strand,
}

impl Step {
    fn new(segment: impl Into<String>, strand: Strand) -> Self {
        Self {
            segment: segment.into(),
            strand,
        }
    }

    fn reverse(&self) -> Self {
        Self {
            segment: self.segment.clone(),
            strand: self.strand.flip(),
        }
    }

    pub fn token(&self) -> String {
        format!("{}{}", self.strand.token_prefix(), self.segment)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SegmentRecord {
    pub name: String,
    pub sequence: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinkRecord {
    pub source: Step,
    pub target: Step,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathRecord {
    pub name: String,
    pub sample: String,
    pub steps: Vec<Step>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NativeGfa {
    pub segments: Vec<SegmentRecord>,
    pub links: Vec<LinkRecord>,
    pub paths: Vec<PathRecord>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FlubbleCandidate {
    pub id: u64,
    pub class_id: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FlubbleBoundary {
    pub open_edge: u64,
    pub close_edge: u64,
}

impl FlubbleBoundary {
    fn ordered(left: u64, right: u64) -> Self {
        if left < right {
            Self {
                open_edge: left,
                close_edge: right,
            }
        } else {
            Self {
                open_edge: right,
                close_edge: left,
            }
        }
    }
}

/// Rust port of `PovuLean.Algorithms.Flubble.detectStack`.
///
/// The Lean detector consumes an already-certified candidate stack and cycle
/// class assignment. This function mirrors that semantic boundary exactly: it
/// does not discover the stack or classes, it only applies the canonical
/// close-after-gap rule and deduplicates boundaries while preserving the Lean
/// membership semantics.
pub fn detect_flubble_stack(stack: &[FlubbleCandidate]) -> Vec<FlubbleBoundary> {
    let mut raw = Vec::new();
    for left_idx in 0..stack.len() {
        let left = &stack[left_idx];
        let rest = &stack[left_idx + 1..];
        if rest.is_empty() {
            continue;
        }
        if rest[0].class_id == left.class_id {
            continue;
        }
        if let Some(right) = rest
            .iter()
            .skip(1)
            .find(|candidate| candidate.class_id == left.class_id)
        {
            raw.push(FlubbleBoundary::ordered(left.id, right.id));
        }
    }

    let mut seen = HashSet::new();
    let mut unique = Vec::new();
    for boundary in raw {
        if seen.insert((boundary.open_edge, boundary.close_edge)) {
            unique.push(boundary);
        }
    }
    unique
}

pub fn gfa_to_vcf_document(
    gfa_path: impl AsRef<Path>,
    ref_file: Option<impl AsRef<Path>>,
) -> Result<VcfDocument> {
    let graph = NativeGfa::load(gfa_path.as_ref())?;
    let reference_names = match ref_file {
        Some(path) => read_reference_names(path.as_ref())?,
        None => Vec::new(),
    };
    graph.to_vcf_document(&reference_names)
}

impl NativeGfa {
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let text = fs::read_to_string(path)?;
        Self::parse(&text)
    }

    pub fn parse(text: &str) -> Result<Self> {
        let mut segments = Vec::new();
        let mut links = Vec::new();
        let mut paths = Vec::new();
        let mut segment_names = HashSet::new();
        let mut path_names = HashSet::new();

        for (line_idx, line) in text.lines().enumerate() {
            let line_no = line_idx + 1;
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            match fields.first().copied() {
                Some("H") => {}
                Some("S") => {
                    if fields.len() < 3 {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: S record requires name and sequence"
                        )));
                    }
                    let name = fields[1];
                    let sequence = fields[2];
                    if name.is_empty() {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: segment name must be nonempty"
                        )));
                    }
                    if sequence == "*" {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: unknown segment sequence '*' is unsupported"
                        )));
                    }
                    if !segment_names.insert(name.to_string()) {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: duplicate segment name '{name}'"
                        )));
                    }
                    segments.push(SegmentRecord {
                        name: name.to_string(),
                        sequence: sequence.to_string(),
                    });
                }
                Some("L") => {
                    if fields.len() < 6 {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: L record requires six core fields"
                        )));
                    }
                    if fields[5] != "0M" {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: only 0M GFA links are supported"
                        )));
                    }
                    links.push(LinkRecord {
                        source: Step::new(fields[1], Strand::from_gfa(fields[2])?),
                        target: Step::new(fields[3], Strand::from_gfa(fields[4])?),
                    });
                }
                Some("P") => {
                    if fields.len() < 4 {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: P record requires name, steps, and overlaps"
                        )));
                    }
                    let name = fields[1];
                    if name.is_empty() {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: path name must be nonempty"
                        )));
                    }
                    if !path_names.insert(name.to_string()) {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: duplicate path name '{name}'"
                        )));
                    }
                    let steps = parse_path_steps(fields[2], line_no)?;
                    if steps.is_empty() {
                        return Err(Error::invalid_gfa(format!(
                            "line {line_no}: path must contain at least one step"
                        )));
                    }
                    paths.push(PathRecord {
                        name: name.to_string(),
                        sample: sample_name(name),
                        steps,
                    });
                }
                Some("W") => {
                    return Err(Error::invalid_gfa(format!(
                        "line {line_no}: W-line parsing is not implemented in the native Rust path yet"
                    )));
                }
                Some("C") | Some("J") => {
                    return Err(Error::invalid_gfa(format!(
                        "line {line_no}: GFA {} records are outside the accepted subset",
                        fields[0]
                    )));
                }
                Some(kind) => {
                    return Err(Error::invalid_gfa(format!(
                        "line {line_no}: unsupported GFA record kind '{kind}'"
                    )));
                }
                None => {}
            }
        }

        let graph = Self {
            segments,
            links,
            paths,
        };
        graph.validate()?;
        Ok(graph)
    }

    fn validate(&self) -> Result<()> {
        if self.segments.is_empty() {
            return Err(Error::invalid_gfa("GFA must contain at least one segment"));
        }
        if self.paths.is_empty() {
            return Err(Error::invalid_gfa("GFA must contain at least one path"));
        }

        let segment_names: HashSet<&str> = self.segments.iter().map(|s| s.name.as_str()).collect();
        for link in &self.links {
            if !segment_names.contains(link.source.segment.as_str()) {
                return Err(Error::invalid_gfa(format!(
                    "link source segment '{}' is not declared",
                    link.source.segment
                )));
            }
            if !segment_names.contains(link.target.segment.as_str()) {
                return Err(Error::invalid_gfa(format!(
                    "link target segment '{}' is not declared",
                    link.target.segment
                )));
            }
        }

        let link_set = self.bidirected_link_set();
        for path in &self.paths {
            for step in &path.steps {
                if !segment_names.contains(step.segment.as_str()) {
                    return Err(Error::invalid_gfa(format!(
                        "path '{}' references undeclared segment '{}'",
                        path.name, step.segment
                    )));
                }
            }
            for window in path.steps.windows(2) {
                let edge = (window[0].clone(), window[1].clone());
                if !link_set.contains(&edge) {
                    return Err(Error::invalid_gfa(format!(
                        "path '{}' has no 0M link between {} and {}",
                        path.name,
                        window[0].token(),
                        window[1].token()
                    )));
                }
            }
        }
        Ok(())
    }

    fn bidirected_link_set(&self) -> HashSet<(Step, Step)> {
        let mut links = HashSet::new();
        for link in &self.links {
            links.insert((link.source.clone(), link.target.clone()));
            links.insert((link.target.reverse(), link.source.reverse()));
        }
        links
    }

    pub fn to_vcf_document(&self, reference_names: &[String]) -> Result<VcfDocument> {
        let reference_idx = self.reference_path_index(reference_names)?;
        let reference = &self.paths[reference_idx];
        let sample_paths = self.sample_paths();
        let samples = self.sample_order();

        let calls =
            if let Some(call) = self.try_hairpin_call(reference_idx, &samples, &sample_paths)? {
                vec![call]
            } else {
                self.flubble_calls(reference_idx, &samples, &sample_paths)?
            };

        Ok(VcfDocument::new(samples, calls)
            .with_source("povu-rs-native")
            .with_contigs(vec![Contig::new(
                reference.name.clone(),
                Some(self.path_len(reference)?),
            )]))
    }

    fn reference_path_index(&self, reference_names: &[String]) -> Result<usize> {
        if reference_names.is_empty() {
            return Ok(0);
        }
        self.paths
            .iter()
            .position(|path| {
                reference_names.iter().any(|name| {
                    path.name == *name || path.sample == *name || path.name.starts_with(name)
                })
            })
            .ok_or_else(|| {
                Error::invalid_gfa(format!(
                    "none of the requested reference names matched a path: {}",
                    reference_names.join(",")
                ))
            })
    }

    fn sample_paths(&self) -> BTreeMap<String, Vec<usize>> {
        let mut by_sample = BTreeMap::<String, Vec<usize>>::new();
        for (idx, path) in self.paths.iter().enumerate() {
            by_sample.entry(path.sample.clone()).or_default().push(idx);
        }
        by_sample
    }

    fn sample_order(&self) -> Vec<String> {
        let mut seen = HashSet::new();
        let mut samples = Vec::new();
        for path in &self.paths {
            if seen.insert(path.sample.as_str()) {
                samples.push(path.sample.clone());
            }
        }
        samples
    }

    fn flubble_calls(
        &self,
        reference_idx: usize,
        samples: &[String],
        sample_paths: &BTreeMap<String, Vec<usize>>,
    ) -> Result<Vec<VariantCall>> {
        let reference = &self.paths[reference_idx];
        let reference_positions = self.reference_positions(reference)?;
        let candidates = self.candidate_sites(reference_idx);
        let leaf_keys = leaf_site_keys(&candidates);
        let mut calls = Vec::new();

        for key in leaf_keys {
            let site = candidates.get(&key).expect("leaf key is from candidates");
            let ref_internal = reference.steps[site.start_idx + 1..site.end_idx].to_vec();
            let start = &reference.steps[site.start_idx];
            let end = &reference.steps[site.end_idx];
            let boundary_id = format!("{}{}", start.token(), end.token());
            let parent_id = nearest_parent_id(&candidates, site, reference);
            let level = nesting_level(&candidates, site) as u32;

            let mut alt_steps = BTreeMap::<Vec<Step>, (usize, u64)>::new();
            let mut sample_alleles = BTreeMap::<String, Vec<GenotypeAllele>>::new();

            for sample in samples {
                let mut alleles = Vec::new();
                for path_idx in &sample_paths[sample] {
                    let path = &self.paths[*path_idx];
                    match allele_steps(path, start, end) {
                        None => alleles.push(GenotypeAllele::Missing),
                        Some(steps) if steps == ref_internal => alleles.push(GenotypeAllele::Ref),
                        Some(steps) => {
                            let next_index = alt_steps.len() + 1;
                            let entry = alt_steps.entry(steps).or_insert((next_index, 0));
                            entry.1 += 1;
                            alleles.push(GenotypeAllele::Alt(entry.0));
                        }
                    }
                }
                sample_alleles.insert(sample.clone(), alleles);
            }

            if alt_steps.is_empty() {
                continue;
            }

            let mut alternates = Vec::new();
            let mut variant_type = None;
            let mut alt_entries: Vec<_> = alt_steps.iter().collect();
            alt_entries.sort_by_key(|(_steps, (index, _count))| *index);
            for (steps, (_index, count)) in alt_entries {
                let alt = self.alternate_for_site(start, &ref_internal, steps, *count)?;
                let alt_type = alt.variant_type();
                if let Some(existing) = variant_type {
                    if existing != alt_type {
                        return Err(Error::unsupported(
                            "native Rust VCF extraction does not yet emit mixed variant types at one site",
                        ));
                    }
                } else {
                    variant_type = Some(alt_type);
                }
                alternates.push(alt);
            }

            let variant_type = variant_type.expect("alt_steps is nonempty");
            let ref_allele = self.ref_allele_for_site(start, &ref_internal, variant_type)?;
            let ref_traversal = ref_traversal_for_site(start, &ref_internal, variant_type);
            let pos = site_position(
                site.start_idx,
                &ref_internal,
                &reference_positions,
                variant_type,
            );
            let genotypes = samples
                .iter()
                .map(|sample| {
                    GenotypeColumn::new(
                        sample.clone(),
                        sample_alleles
                            .remove(sample)
                            .expect("genotype alleles were populated for each sample"),
                    )
                })
                .collect::<Vec<_>>();
            let reference_allele_count = genotypes
                .iter()
                .flat_map(|column| column.alleles.iter())
                .filter(|allele| matches!(allele, GenotypeAllele::Ref))
                .count() as u64;

            let source = match parent_id {
                Some(parent) => VariantSource::nested_flubble(
                    boundary_id.clone(),
                    parent,
                    level,
                    step_order_key(start, site.start_idx),
                    step_order_key(end, site.end_idx),
                ),
                None => VariantSource::flubble(
                    boundary_id.clone(),
                    level,
                    step_order_key(start, site.start_idx),
                    step_order_key(end, site.end_idx),
                ),
            };

            calls.push(VariantCall {
                source,
                chrom: reference.name.clone(),
                contig_order: 0,
                pos,
                id: boundary_id,
                ref_allele,
                ref_traversal,
                alternates,
                variant_type,
                tangled: false,
                reference_allele_count,
                genotypes,
            });
        }

        calls.sort_by_key(VariantCall::order_key);
        Ok(calls)
    }

    fn candidate_sites(&self, reference_idx: usize) -> BTreeMap<SiteKey, SiteDraft> {
        let reference = &self.paths[reference_idx];
        let mut candidates = BTreeMap::<SiteKey, SiteDraft>::new();
        for (path_idx, path) in self.paths.iter().enumerate() {
            if path_idx == reference_idx {
                continue;
            }
            let matches = collinear_matches(&reference.steps, &path.steps);
            for window in matches.windows(2) {
                let (ref_start, path_start) = window[0];
                let (ref_end, path_end) = window[1];
                if ref_end <= ref_start || path_end <= path_start {
                    continue;
                }
                let ref_internal = &reference.steps[ref_start + 1..ref_end];
                let alt_internal = &path.steps[path_start + 1..path_end];
                if ref_internal == alt_internal {
                    continue;
                }
                candidates.entry((ref_start, ref_end)).or_insert(SiteDraft {
                    start_idx: ref_start,
                    end_idx: ref_end,
                });
            }
        }
        candidates
    }

    fn try_hairpin_call(
        &self,
        reference_idx: usize,
        samples: &[String],
        sample_paths: &BTreeMap<String, Vec<usize>>,
    ) -> Result<Option<VariantCall>> {
        let reference = &self.paths[reference_idx];
        let Some((alt_idx, _)) = self
            .paths
            .iter()
            .enumerate()
            .find(|(idx, path)| *idx != reference_idx && is_reverse_path(reference, path))
        else {
            return Ok(None);
        };

        let positions = self.reference_positions(reference)?;
        let Some(first) = reference.steps.first() else {
            return Ok(None);
        };
        let Some(last) = reference.steps.last() else {
            return Ok(None);
        };
        let boundary_id = format!("{}{}", first.token(), last.token());
        let ref_allele = self.steps_sequence(&reference.steps)?;
        let alt_allele = self.steps_sequence(&self.paths[alt_idx].steps)?;
        let ref_traversal = traversal(&reference.steps);
        let alt_traversal = traversal(&self.paths[alt_idx].steps);

        let genotypes = samples
            .iter()
            .map(|sample| {
                let alleles = sample_paths[sample]
                    .iter()
                    .map(|path_idx| {
                        if *path_idx == reference_idx {
                            GenotypeAllele::Ref
                        } else if *path_idx == alt_idx {
                            GenotypeAllele::Alt(1)
                        } else {
                            GenotypeAllele::Missing
                        }
                    })
                    .collect();
                GenotypeColumn::new(sample.clone(), alleles)
            })
            .collect::<Vec<_>>();
        let reference_allele_count = genotypes
            .iter()
            .flat_map(|column| column.alleles.iter())
            .filter(|allele| matches!(allele, GenotypeAllele::Ref))
            .count() as u64;

        Ok(Some(VariantCall {
            source: VariantSource::hairpin(
                boundary_id.clone(),
                step_order_key(first, 0),
                step_order_key(last, reference.steps.len().saturating_sub(1)),
            ),
            chrom: reference.name.clone(),
            contig_order: 0,
            pos: positions.get(1).copied().unwrap_or(1),
            id: boundary_id,
            ref_allele: ref_allele.clone(),
            ref_traversal,
            alternates: vec![AlternateAllele::new(
                AlleleConstruction::reverse_substitution(ref_allele, alt_allele),
                alt_traversal,
                1,
            )],
            variant_type: VariantType::Subr,
            tangled: false,
            reference_allele_count,
            genotypes,
        }))
    }

    fn alternate_for_site(
        &self,
        start: &Step,
        ref_internal: &[Step],
        alt_internal: &[Step],
        count: u64,
    ) -> Result<AlternateAllele> {
        let anchor = self.step_sequence(start)?;
        let alt_traversal;
        let construction = if ref_internal.is_empty() {
            let inserted = self.steps_sequence(alt_internal)?;
            alt_traversal = traversal_with_anchor(start, alt_internal);
            AlleleConstruction::insertion(anchor, inserted)
        } else if alt_internal.is_empty() {
            let deleted = self.steps_sequence(ref_internal)?;
            alt_traversal = start.token();
            AlleleConstruction::deletion(anchor, deleted)
        } else {
            let reference = self.steps_sequence(ref_internal)?;
            let alternate = self.steps_sequence(alt_internal)?;
            alt_traversal = traversal(alt_internal);
            AlleleConstruction::substitution(reference, alternate)
        };
        Ok(AlternateAllele::new(construction, alt_traversal, count))
    }

    fn ref_allele_for_site(
        &self,
        start: &Step,
        ref_internal: &[Step],
        variant_type: VariantType,
    ) -> Result<String> {
        match variant_type {
            VariantType::Ins => self.step_sequence(start),
            VariantType::Del => Ok(format!(
                "{}{}",
                self.step_sequence(start)?,
                self.steps_sequence(ref_internal)?
            )),
            VariantType::Sub => self.steps_sequence(ref_internal),
            VariantType::Subr => unreachable!("SUBR is handled by hairpin extraction"),
        }
    }

    fn step_sequence(&self, step: &Step) -> Result<String> {
        let segment = self
            .segments
            .iter()
            .find(|segment| segment.name == step.segment)
            .ok_or_else(|| Error::invalid_gfa(format!("unknown segment '{}'", step.segment)))?;
        Ok(match step.strand {
            Strand::Forward => segment.sequence.clone(),
            Strand::Reverse => reverse_complement(&segment.sequence),
        })
    }

    fn steps_sequence(&self, steps: &[Step]) -> Result<String> {
        let mut sequence = String::new();
        for step in steps {
            sequence.push_str(&self.step_sequence(step)?);
        }
        Ok(sequence)
    }

    fn path_len(&self, path: &PathRecord) -> Result<u64> {
        Ok(self.steps_sequence(&path.steps)?.len() as u64)
    }

    fn reference_positions(&self, reference: &PathRecord) -> Result<Vec<u64>> {
        let mut positions = Vec::with_capacity(reference.steps.len());
        let mut next = 1_u64;
        for step in &reference.steps {
            positions.push(next);
            next += self.step_sequence(step)?.len() as u64;
        }
        Ok(positions)
    }
}

type SiteKey = (usize, usize);

#[derive(Clone, Debug, PartialEq, Eq)]
struct SiteDraft {
    start_idx: usize,
    end_idx: usize,
}

fn parse_path_steps(value: &str, line_no: usize) -> Result<Vec<Step>> {
    if value.is_empty() {
        return Ok(Vec::new());
    }
    value
        .split(',')
        .map(|token| {
            let Some(strand_char) = token.chars().last() else {
                return Err(Error::invalid_gfa(format!(
                    "line {line_no}: empty path step"
                )));
            };
            let segment = &token[..token.len() - strand_char.len_utf8()];
            if segment.is_empty() {
                return Err(Error::invalid_gfa(format!(
                    "line {line_no}: path step has empty segment name"
                )));
            }
            Ok(Step::new(segment, Strand::from_step_suffix(strand_char)?))
        })
        .collect()
}

fn sample_name(path_name: &str) -> String {
    path_name
        .split('#')
        .next()
        .filter(|name| !name.is_empty())
        .unwrap_or(path_name)
        .to_string()
}

fn read_reference_names(path: &Path) -> Result<Vec<String>> {
    let text = fs::read_to_string(path)?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| line.split_whitespace().next().unwrap_or(line).to_string())
        .collect())
}

fn collinear_matches(reference: &[Step], path: &[Step]) -> Vec<(usize, usize)> {
    let rows = reference.len();
    let cols = path.len();
    let mut dp = vec![vec![0_usize; cols + 1]; rows + 1];
    for i in (0..rows).rev() {
        for j in (0..cols).rev() {
            dp[i][j] = if reference[i] == path[j] {
                1 + dp[i + 1][j + 1]
            } else {
                dp[i + 1][j].max(dp[i][j + 1])
            };
        }
    }

    let mut matches = Vec::new();
    let mut i = 0;
    let mut j = 0;
    while i < rows && j < cols {
        if reference[i] == path[j] {
            matches.push((i, j));
            i += 1;
            j += 1;
        } else if dp[i + 1][j] >= dp[i][j + 1] {
            i += 1;
        } else {
            j += 1;
        }
    }
    matches
}

fn leaf_site_keys(candidates: &BTreeMap<SiteKey, SiteDraft>) -> Vec<SiteKey> {
    candidates
        .iter()
        .filter_map(|(key, site)| {
            let contains_child = candidates
                .values()
                .any(|other| site.start_idx < other.start_idx && other.end_idx < site.end_idx);
            (!contains_child).then_some(*key)
        })
        .collect()
}

fn nesting_level(candidates: &BTreeMap<SiteKey, SiteDraft>, site: &SiteDraft) -> usize {
    candidates
        .values()
        .filter(|other| other.start_idx < site.start_idx && site.end_idx < other.end_idx)
        .count()
}

fn nearest_parent_id(
    candidates: &BTreeMap<SiteKey, SiteDraft>,
    site: &SiteDraft,
    reference: &PathRecord,
) -> Option<String> {
    candidates
        .values()
        .filter(|other| other.start_idx < site.start_idx && site.end_idx < other.end_idx)
        .min_by_key(|other| other.end_idx - other.start_idx)
        .map(|parent| parent_boundary_id(parent, reference))
}

fn parent_boundary_id(site: &SiteDraft, reference: &PathRecord) -> String {
    format!(
        "{}{}",
        reference.steps[site.start_idx].token(),
        reference.steps[site.end_idx].token()
    )
}

fn allele_steps(path: &PathRecord, start: &Step, end: &Step) -> Option<Vec<Step>> {
    let start_idx = path.steps.iter().position(|step| step == start)?;
    let end_idx = path.steps[start_idx + 1..]
        .iter()
        .position(|step| step == end)
        .map(|offset| start_idx + 1 + offset)?;
    Some(path.steps[start_idx + 1..end_idx].to_vec())
}

fn ref_traversal_for_site(
    start: &Step,
    ref_internal: &[Step],
    variant_type: VariantType,
) -> String {
    match variant_type {
        VariantType::Ins => start.token(),
        VariantType::Del => traversal_with_anchor(start, ref_internal),
        VariantType::Sub => traversal(ref_internal),
        VariantType::Subr => unreachable!("SUBR is handled by hairpin extraction"),
    }
}

fn site_position(
    start_idx: usize,
    ref_internal: &[Step],
    reference_positions: &[u64],
    variant_type: VariantType,
) -> u64 {
    match variant_type {
        VariantType::Ins | VariantType::Del => reference_positions[start_idx],
        VariantType::Sub => {
            if ref_internal.is_empty() {
                reference_positions[start_idx]
            } else {
                reference_positions[start_idx + 1]
            }
        }
        VariantType::Subr => unreachable!("SUBR is handled by hairpin extraction"),
    }
}

fn traversal(steps: &[Step]) -> String {
    steps.iter().map(Step::token).collect()
}

fn traversal_with_anchor(start: &Step, rest: &[Step]) -> String {
    let mut value = start.token();
    value.push_str(&traversal(rest));
    value
}

fn is_reverse_path(reference: &PathRecord, path: &PathRecord) -> bool {
    reference.steps.len() == path.steps.len()
        && reference
            .steps
            .iter()
            .zip(path.steps.iter().rev())
            .all(|(left, right)| {
                left.segment == right.segment && left.strand == right.strand.flip()
            })
}

fn step_order_key(step: &Step, fallback: usize) -> u64 {
    step.segment.parse::<u64>().unwrap_or(fallback as u64)
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' | 'a' => 'T',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'T' | 't' => 'A',
            'N' | 'n' => 'N',
            other => other,
        })
        .collect()
}
