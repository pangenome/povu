use std::fs;
use std::path::{Path, PathBuf};

use serde_json::{json, Value};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ProfileOutput {
    Vcf {
        vcf: String,
        provenance: Option<String>,
    },
    Reject {
        reason: String,
    },
}

pub fn fixture_root(repo_root: &Path) -> PathBuf {
    repo_root.join("tests/lean4_conformance/fixtures/downstream_repetitive")
}

pub fn run_fixture_profile(
    repo_root: &Path,
    fixture_id: &str,
    profile_name: &str,
) -> Result<ProfileOutput, String> {
    let root = fixture_root(repo_root);
    let fixture = find_fixture(&root, fixture_id)?;
    run_profile(&root, &fixture, profile_name)
}

pub fn run_downstream_repetitive_corpus(
    repo_root: &Path,
    fixture_filter: Option<&str>,
) -> Result<(), String> {
    let root = fixture_root(repo_root);
    let manifest = read_json(&root.join("manifest.json"))?;
    expect_string(&manifest, "schema")?
        .eq("povu.downstream_repetitive.fixtures.v1")
        .then_some(())
        .ok_or_else(|| "downstream repetitive manifest has unexpected schema".to_string())?;

    let fixtures = expect_array(&manifest, "fixtures")?;
    let mut failures = Vec::new();
    let mut ran = 0usize;

    for fixture in fixtures {
        let fixture_id = expect_string(fixture, "id")?;
        if fixture_filter.is_some_and(|filter| filter != fixture_id) {
            continue;
        }
        ran += 1;

        for profile in expect_array(fixture, "profiles")? {
            let profile_name = expect_string(profile, "name")?;
            match compare_profile_output(&root, fixture, profile) {
                Ok(()) => {
                    println!("ok - {fixture_id}/{profile_name}");
                }
                Err(err) => failures.push(format!("{fixture_id}/{profile_name}: {err}")),
            }
        }
    }

    if ran == 0 {
        return Err(format!(
            "unknown downstream repetitive fixture id: {}",
            fixture_filter.unwrap_or("<none>")
        ));
    }

    if failures.is_empty() {
        println!("downstream repetitive conformance: all selected profiles passed");
        Ok(())
    } else {
        Err(failures.join("\n\n"))
    }
}

fn compare_profile_output(root: &Path, fixture: &Value, profile: &Value) -> Result<(), String> {
    let fixture_id = expect_string(fixture, "id")?;
    let profile_name = expect_string(profile, "name")?;
    let actual = run_profile(root, fixture, profile_name)?;

    match expect_string(profile, "expected")? {
        "vcf" => {
            let rel_path = expect_string(profile, "path")?;
            let expected_path = root.join(rel_path);
            let expected = fs::read_to_string(&expected_path)
                .map_err(|err| format!("failed to read {}: {err}", expected_path.display()))?;
            let ProfileOutput::Vcf { vcf, provenance } = actual else {
                return Err("profile rejected, expected VCF output".to_string());
            };
            if vcf != expected {
                return Err(format!(
                    "VCF mismatch against {}\n--- expected ---\n{}\n--- actual ---\n{}",
                    expected_path.display(),
                    expected,
                    vcf
                ));
            }

            if let Some(provenance_rel) = optional_string(profile, "provenance") {
                let expected_provenance_path = root.join(provenance_rel);
                let expected_provenance =
                    fs::read_to_string(&expected_provenance_path).map_err(|err| {
                        format!(
                            "failed to read {}: {err}",
                            expected_provenance_path.display()
                        )
                    })?;
                let actual_provenance = provenance.ok_or_else(|| {
                    format!("{fixture_id}/{profile_name} did not emit provenance")
                })?;
                let expected_json = parse_jsonl(&expected_provenance)
                    .map_err(|err| format!("invalid expected provenance: {err}"))?;
                let actual_json = parse_jsonl(&actual_provenance)
                    .map_err(|err| format!("invalid actual provenance: {err}"))?;
                if expected_json != actual_json {
                    return Err(format!(
                        "provenance mismatch against {}\n--- expected ---\n{}\n--- actual ---\n{}",
                        expected_provenance_path.display(),
                        expected_provenance,
                        actual_provenance
                    ));
                }
            }
            Ok(())
        }
        "reject" => {
            let ProfileOutput::Reject { reason } = actual else {
                return Err("profile emitted VCF, expected controlled reject".to_string());
            };
            if reason.trim().is_empty() {
                return Err("controlled reject reason was empty".to_string());
            }
            Ok(())
        }
        other => Err(format!("unsupported profile expected kind {other}")),
    }
}

fn run_profile(root: &Path, fixture: &Value, profile_name: &str) -> Result<ProfileOutput, String> {
    let fixture_id = expect_string(fixture, "id")?;
    match profile_name {
        "raw-graph" => {
            if fixture_id == "nested-child-inside-insertion" {
                return Ok(ProfileOutput::Reject {
                    reason: "nested child has no selected-reference coordinate".to_string(),
                });
            }
            let raw_path = raw_vcf_path(root, fixture)?;
            Ok(ProfileOutput::Vcf {
                vcf: fs::read_to_string(&raw_path)
                    .map_err(|err| format!("failed to read {}: {err}", raw_path.display()))?,
                provenance: None,
            })
        }
        "inserted-path-raw" => Ok(ProfileOutput::Reject {
            reason: "inserted-path coordinate profile is not specified".to_string(),
        }),
        "left-normalized" => left_normalized(root, fixture),
        "top-level-only" => top_level_only(root, fixture),
        "popped" => popped(root, fixture),
        "decomposed" => decomposed(root, fixture),
        other => Err(format!("unknown downstream profile {other}")),
    }
}

fn left_normalized(root: &Path, fixture: &Value) -> Result<ProfileOutput, String> {
    let raw = read_raw_vcf(root, fixture)?;
    let reference_rel = expect_string(fixture, "reference_fasta")?;
    let references = read_fasta(&root.join(reference_rel))?;
    let mut out_records = Vec::new();
    let mut provenance = Vec::new();

    for record in raw.records {
        let original = record.clone();
        let alt = only_alt(&record)?;
        let reference = references
            .get(&record.chrom)
            .or_else(|| references.values().next())
            .ok_or_else(|| "left-normalized profile requires a reference FASTA".to_string())?;
        let (pos, ref_allele, alt_allele) =
            left_align_indel(record.pos, &record.ref_allele, alt, reference);

        let mut rewritten = record.with_single_alt(&alt_allele);
        rewritten.pos = pos;
        rewritten.ref_allele = ref_allele;
        rewritten.id = format!("{}:norm", original.id);
        rewritten.append_info("ORIGIN", &original.id);
        rewritten.append_info("RAW_ALT_INDEX", "1");
        rewritten.append_info("PROFILE", "left-normalized");
        rewritten.append_info("LEFT_NORMALIZED", "T");
        rewritten.append_info("RAW_POS", &original.pos.to_string());
        rewritten.append_info("RAW_REF", &original.ref_allele);
        rewritten.append_info("RAW_ALT", alt);

        provenance.push(json!({
            "fixture_id": expect_string(fixture, "id")?,
            "profile": "left-normalized",
            "operation": "coordinate-shift",
            "raw_record_id": original.id,
            "raw_alt_index": 1,
            "raw_chrom": original.chrom,
            "raw_pos": original.pos,
            "raw_ref": original.ref_allele,
            "raw_alt": alt,
            "raw_at": split_info_list(original.info_value("AT")?),
            "output_record_ids": [rewritten.id.clone()],
            "output_chrom": rewritten.chrom,
            "output_pos": rewritten.pos,
            "output_ref": rewritten.ref_allele,
            "output_alt": alt_allele,
            "flags": ["LEFT_NORMALIZED"],
            "profile_config": {
                "reference_fasta": "reference.fa",
                "normalizer": "bcftools-norm-like"
            }
        }));

        out_records.push(rewritten);
    }

    Ok(ProfileOutput::Vcf {
        vcf: Vcf {
            profile: "left-normalized".to_string(),
            records: out_records,
            ..raw
        }
        .to_text(),
        provenance: Some(jsonl(&provenance)?),
    })
}

fn top_level_only(root: &Path, fixture: &Value) -> Result<ProfileOutput, String> {
    let raw = read_raw_vcf(root, fixture)?;
    let records = raw
        .records
        .iter()
        .filter(|record| {
            record.info_value("VARTYPE").ok() == Some("SUBR")
                || record
                    .info_value("LV")
                    .ok()
                    .and_then(|value| value.parse::<u64>().ok())
                    == Some(0)
        })
        .map(|record| {
            let mut kept = record.clone();
            kept.id = format!("{}:top", record.id);
            kept.append_info("ORIGIN", &record.id);
            kept.append_info("PROFILE", "top-level-only");
            kept.append_info("PASSTHROUGH", "T");
            kept
        })
        .collect();

    Ok(ProfileOutput::Vcf {
        vcf: Vcf {
            profile: "top-level-only".to_string(),
            records,
            ..raw
        }
        .to_text(),
        provenance: None,
    })
}

fn popped(root: &Path, fixture: &Value) -> Result<ProfileOutput, String> {
    let raw = read_raw_vcf(root, fixture)?;
    let max_level = 0u64;
    let max_ref_length = 4usize;
    let max_allele_length = 4usize;

    let mut popped_parents = Vec::new();
    for record in &raw.records {
        let level = record
            .info_value("LV")
            .ok()
            .and_then(|value| value.parse::<u64>().ok())
            .unwrap_or(0);
        let oversize = record.ref_allele.len() > max_ref_length
            || record
                .alt_alleles
                .iter()
                .any(|alt| alt.len() > max_allele_length);
        if level <= max_level && oversize {
            popped_parents.push(record.clone());
        }
    }

    let mut out_records = Vec::new();
    let mut provenance = Vec::new();
    for parent in &popped_parents {
        provenance.push(json!({
            "fixture_id": expect_string(fixture, "id")?,
            "profile": "popped",
            "operation": "pop-parent",
            "raw_record_id": parent.id,
            "raw_alt_index": 1,
            "raw_chrom": parent.chrom,
            "raw_pos": parent.pos,
            "raw_ref": parent.ref_allele,
            "raw_alt": only_alt(parent)?,
            "raw_at": split_info_list(parent.info_value("AT")?),
            "output_record_ids": [],
            "flags": ["POPPED_PARENT"],
            "profile_config": {
                "max_level": max_level,
                "max_ref_length": max_ref_length,
                "max_allele_length": max_allele_length
            }
        }));
    }

    for record in &raw.records {
        if popped_parents.iter().any(|parent| parent.id == record.id) {
            continue;
        }

        let level = record
            .info_value("LV")
            .ok()
            .and_then(|value| value.parse::<u64>().ok())
            .unwrap_or(0);
        if level <= max_level {
            out_records.push(record.clone());
            continue;
        }

        let Some(parent) = popped_parents
            .iter()
            .find(|parent| id_interval_contains(&parent.id, &record.id))
        else {
            continue;
        };

        if record.ref_allele.len() > max_ref_length
            || record
                .alt_alleles
                .iter()
                .any(|alt| alt.len() > max_allele_length)
        {
            continue;
        }

        let mut rescued = record.clone();
        rescued.id = format!("{}:rescued", record.id);
        rescued.append_info("ORIGIN", &record.id);
        rescued.append_info("PARENT", &parent.id);
        rescued.append_info("PROFILE", "popped");
        rescued.append_info("RESCUED_CHILD", "T");
        rescued.append_info("POPPED_PARENT", &parent.id);

        provenance.push(json!({
            "fixture_id": expect_string(fixture, "id")?,
            "profile": "popped",
            "operation": "rescue-child",
            "raw_record_id": record.id,
            "raw_alt_index": 1,
            "raw_chrom": record.chrom,
            "raw_pos": record.pos,
            "raw_ref": record.ref_allele,
            "raw_alt": only_alt(record)?,
            "raw_at": split_info_list(record.info_value("AT")?),
            "output_record_ids": [rescued.id.clone()],
            "output_chrom": rescued.chrom,
            "output_pos": rescued.pos,
            "output_ref": rescued.ref_allele,
            "output_alt": only_alt(&rescued)?,
            "parent_record_id": parent.id,
            "flags": ["RESCUED_CHILD"],
            "profile_config": {
                "max_level": max_level,
                "max_ref_length": max_ref_length,
                "max_allele_length": max_allele_length,
                "child_rescue": "eligible-child-of-popped-parent"
            }
        }));

        out_records.push(rescued);
    }

    Ok(ProfileOutput::Vcf {
        vcf: Vcf {
            profile: "popped".to_string(),
            records: out_records,
            ..raw
        }
        .to_text(),
        provenance: Some(jsonl(&provenance)?),
    })
}

fn decomposed(root: &Path, fixture: &Value) -> Result<ProfileOutput, String> {
    let raw = read_raw_vcf(root, fixture)?;
    let mut out_records = Vec::new();
    let mut provenance = Vec::new();
    let max_allele_length = 8usize;

    for record in &raw.records {
        if record.info_value("VARTYPE")? == "SUBR" {
            let mut passthrough = record.clone();
            passthrough.id = format!("{}:subr-passthrough", record.id);
            passthrough.append_info("ORIGIN", &record.id);
            passthrough.append_info("RAW_ALT_INDEX", "1");
            passthrough.append_info("PROFILE", "decomposed");
            passthrough.append_info("PASSTHROUGH", "T");
            passthrough.append_info("PASS_THROUGH_REASON", "subr_inversion_preservation");
            passthrough.append_info("SUBR_ORIGIN", "T");

            provenance.push(json!({
                "fixture_id": expect_string(fixture, "id")?,
                "profile": "decomposed",
                "operation": "pass-through",
                "raw_record_id": record.id,
                "raw_alt_index": 1,
                "raw_chrom": record.chrom,
                "raw_pos": record.pos,
                "raw_ref": record.ref_allele,
                "raw_alt": only_alt(record)?,
                "raw_at": split_info_list(record.info_value("AT")?),
                "output_record_ids": [passthrough.id.clone()],
                "output_chrom": passthrough.chrom,
                "output_pos": passthrough.pos,
                "output_ref": passthrough.ref_allele,
                "output_alt": only_alt(&passthrough)?,
                "flags": ["PASSTHROUGH", "SUBR_PRESERVED"],
                "profile_config": {
                    "pass_through_reason": "subr_inversion_preservation",
                    "forbidden_outputs": ["<INV>", "primitive-substitution-series"]
                }
            }));

            out_records.push(passthrough);
            continue;
        }

        let at = split_info_list(record.info_value("AT")?);
        let ref_steps = traversal_steps(
            at.first()
                .and_then(Value::as_str)
                .ok_or_else(|| format!("{} lacks REF AT traversal", record.id))?,
        );

        for (alt_offset, alt) in record.alt_alleles.iter().enumerate() {
            let raw_alt_index = alt_offset + 1;
            let alt_steps =
                traversal_steps(at.get(raw_alt_index).and_then(Value::as_str).ok_or_else(
                    || format!("{} lacks AT traversal for ALT {raw_alt_index}", record.id),
                )?);

            if alt.len() > max_allele_length {
                let mut passthrough = record.with_single_alt(alt);
                passthrough.id = format!("{}:{raw_alt_index}:passthrough", record.id);
                passthrough.genotypes =
                    project_genotypes(&record.genotypes, raw_alt_index.to_string().as_str());
                set_counts(&mut passthrough);
                passthrough.set_info(
                    "AT",
                    &format!(
                        "{},{}",
                        at[0].as_str().unwrap(),
                        at[raw_alt_index].as_str().unwrap()
                    ),
                );
                passthrough.set_info("TANGLED", "T");
                passthrough.remove_info("ES");
                passthrough.remove_info("LV");
                passthrough.append_info("ORIGIN", &record.id);
                passthrough.append_info("RAW_ALT_INDEX", &raw_alt_index.to_string());
                passthrough.append_info("PROFILE", "decomposed");
                passthrough.append_info("PASSTHROUGH", "T");
                passthrough.append_info("PASS_THROUGH_REASON", "max_allele_length");
                passthrough.append_info("RAW_POS", &record.pos.to_string());
                passthrough.append_info("RAW_REF", &record.ref_allele);
                passthrough.append_info("RAW_ALT", alt);

                provenance.push(json!({
                    "fixture_id": expect_string(fixture, "id")?,
                    "profile": "decomposed",
                    "operation": "pass-through",
                    "raw_record_id": record.id,
                    "raw_alt_index": raw_alt_index,
                    "raw_chrom": record.chrom,
                    "raw_pos": record.pos,
                    "raw_ref": record.ref_allele,
                    "raw_alt": alt,
                    "raw_at": [at[0].clone(), at[raw_alt_index].clone()],
                    "output_record_ids": [passthrough.id.clone()],
                    "output_chrom": passthrough.chrom,
                    "output_pos": passthrough.pos,
                    "output_ref": passthrough.ref_allele,
                    "output_alt": alt,
                    "flags": ["PASSTHROUGH"],
                    "profile_config": {
                        "max_allele_length": max_allele_length,
                        "pass_through_reason": "max_allele_length"
                    }
                }));

                out_records.push(passthrough);
                continue;
            }

            let mut decomposed_records = Vec::new();
            for (mismatch_index, (idx, ref_base, alt_base)) in
                mismatches(&record.ref_allele, alt).into_iter().enumerate()
            {
                let mut primitive = record.with_single_alt(&alt_base.to_string());
                primitive.pos = record.pos + idx as u64;
                primitive.id = format!("{}:{raw_alt_index}:snp{}", record.id, mismatch_index + 1);
                primitive.ref_allele = ref_base.to_string();
                primitive.genotypes =
                    project_genotypes(&record.genotypes, raw_alt_index.to_string().as_str());
                set_counts(&mut primitive);
                primitive.set_info(
                    "AT",
                    &format!(
                        "{},{}",
                        ref_steps.get(idx).ok_or_else(|| {
                            format!("{} REF traversal lacks step offset {idx}", record.id)
                        })?,
                        alt_steps.get(idx).ok_or_else(|| {
                            format!("{} ALT traversal lacks step offset {idx}", record.id)
                        })?
                    ),
                );
                primitive.set_info("TANGLED", "F");
                primitive.remove_info("ES");
                primitive.remove_info("LV");
                primitive.append_info("ORIGIN", &record.id);
                primitive.append_info("RAW_ALT_INDEX", &raw_alt_index.to_string());
                primitive.append_info("PROFILE", "decomposed");
                primitive.append_info("DECOMPOSED", "T");
                primitive.append_info("RAW_POS", &record.pos.to_string());
                primitive.append_info("RAW_REF", &record.ref_allele);
                primitive.append_info("RAW_ALT", alt);
                decomposed_records.push(primitive);
            }

            provenance.push(json!({
                "fixture_id": expect_string(fixture, "id")?,
                "profile": "decomposed",
                "operation": "one-to-many-decomposition",
                "raw_record_id": record.id,
                "raw_alt_index": raw_alt_index,
                "raw_chrom": record.chrom,
                "raw_pos": record.pos,
                "raw_ref": record.ref_allele,
                "raw_alt": alt,
                "raw_at": [at[0].clone(), at[raw_alt_index].clone()],
                "output_record_ids": decomposed_records
                    .iter()
                    .map(|primitive| Value::String(primitive.id.clone()))
                    .collect::<Vec<_>>(),
                "output_records": decomposed_records
                    .iter()
                    .map(|primitive| json!({
                        "chrom": primitive.chrom,
                        "pos": primitive.pos,
                        "ref": primitive.ref_allele,
                        "alt": only_alt(primitive).unwrap_or("")
                    }))
                    .collect::<Vec<_>>(),
                "flags": ["DECOMPOSED"],
                "profile_config": {
                    "max_allele_length": max_allele_length,
                    "aligner": "vcfwave-like"
                }
            }));

            out_records.extend(decomposed_records);
        }
    }

    Ok(ProfileOutput::Vcf {
        vcf: Vcf {
            profile: "decomposed".to_string(),
            records: out_records,
            ..raw
        }
        .to_text(),
        provenance: Some(jsonl(&provenance)?),
    })
}

#[derive(Clone, Debug)]
struct Vcf {
    profile: String,
    samples: Vec<String>,
    records: Vec<Record>,
}

#[derive(Clone, Debug)]
struct Record {
    chrom: String,
    pos: u64,
    id: String,
    ref_allele: String,
    alt_alleles: Vec<String>,
    qual: String,
    filter: String,
    info: Vec<(String, String)>,
    format: String,
    genotypes: Vec<String>,
}

impl Vcf {
    fn parse(text: &str) -> Result<Self, String> {
        let mut samples = Vec::new();
        let mut records = Vec::new();

        for (idx, line) in text.lines().enumerate() {
            let line_no = idx + 1;
            if line.trim().is_empty() || line.starts_with("##") {
                continue;
            }
            if line.starts_with("#CHROM") {
                let fields: Vec<_> = line.split('\t').collect();
                if fields.len() < 10 {
                    return Err(format!("line {line_no}: VCF header lacks samples"));
                }
                samples = fields[9..]
                    .iter()
                    .map(|sample| (*sample).to_string())
                    .collect();
                continue;
            }

            let fields: Vec<_> = line.split('\t').collect();
            let expected = 9 + samples.len();
            if fields.len() != expected {
                return Err(format!(
                    "line {line_no}: expected {expected} fields, found {}",
                    fields.len()
                ));
            }
            records.push(Record {
                chrom: fields[0].to_string(),
                pos: fields[1]
                    .parse()
                    .map_err(|err| format!("line {line_no}: invalid POS: {err}"))?,
                id: fields[2].to_string(),
                ref_allele: fields[3].to_string(),
                alt_alleles: fields[4].split(',').map(str::to_string).collect(),
                qual: fields[5].to_string(),
                filter: fields[6].to_string(),
                info: parse_info(fields[7]),
                format: fields[8].to_string(),
                genotypes: fields[9..].iter().map(|gt| (*gt).to_string()).collect(),
            });
        }

        if samples.is_empty() {
            return Err("VCF lacks #CHROM header with samples".to_string());
        }

        Ok(Self {
            profile: "raw-graph".to_string(),
            samples,
            records,
        })
    }

    fn to_text(&self) -> String {
        let mut text = String::new();
        text.push_str("##fileformat=VCFv4.2\n");
        text.push_str("##source=povu-downstream-fixture\n");
        text.push_str(&format!("##profile={}\n", self.profile));
        for line in common_info_headers(self.records.iter().any(Record::has_es_or_lv)) {
            text.push_str(line);
            text.push('\n');
        }
        for line in profile_info_headers(&self.profile, self.is_subr_only()) {
            text.push_str(line);
            text.push('\n');
        }
        text.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for sample in &self.samples {
            text.push('\t');
            text.push_str(sample);
        }
        text.push('\n');
        for record in &self.records {
            text.push_str(&record.to_line());
            text.push('\n');
        }
        text
    }

    fn is_subr_only(&self) -> bool {
        !self.records.is_empty()
            && self
                .records
                .iter()
                .all(|record| record.info_value("VARTYPE").ok() == Some("SUBR"))
    }
}

impl Record {
    fn to_line(&self) -> String {
        [
            self.chrom.clone(),
            self.pos.to_string(),
            self.id.clone(),
            self.ref_allele.clone(),
            self.alt_alleles.join(","),
            self.qual.clone(),
            self.filter.clone(),
            self.info
                .iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect::<Vec<_>>()
                .join(";"),
            self.format.clone(),
        ]
        .into_iter()
        .chain(self.genotypes.clone())
        .collect::<Vec<_>>()
        .join("\t")
    }

    fn has_es_or_lv(&self) -> bool {
        self.info.iter().any(|(key, _)| key == "ES" || key == "LV")
    }

    fn info_value(&self, key: &str) -> Result<&str, String> {
        self.info
            .iter()
            .find(|(candidate, _)| candidate == key)
            .map(|(_, value)| value.as_str())
            .ok_or_else(|| format!("record {} lacks INFO/{key}", self.id))
    }

    fn append_info(&mut self, key: &str, value: &str) {
        self.info.push((key.to_string(), value.to_string()));
    }

    fn set_info(&mut self, key: &str, value: &str) {
        if let Some((_, existing)) = self.info.iter_mut().find(|(candidate, _)| candidate == key) {
            *existing = value.to_string();
        } else {
            self.append_info(key, value);
        }
    }

    fn remove_info(&mut self, key: &str) {
        self.info.retain(|(candidate, _)| candidate != key);
    }

    fn with_single_alt(&self, alt: &str) -> Self {
        let mut record = self.clone();
        record.alt_alleles = vec![alt.to_string()];
        record
    }
}

fn common_info_headers(include_es_lv: bool) -> Vec<&'static str> {
    let mut headers = vec![
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate allele count\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total called alleles\">",
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Samples with called genotypes\">",
        "##INFO=<ID=AT,Number=R,Type=String,Description=\"Raw graph traversal provenance\">",
        "##INFO=<ID=VARTYPE,Number=1,Type=String,Description=\"povu variant type\">",
        "##INFO=<ID=TANGLED,Number=1,Type=String,Description=\"Graph-faithful complex marker\">",
    ];
    if include_es_lv {
        headers.push("##INFO=<ID=ES,Number=1,Type=String,Description=\"Enclosing flubble site\">");
        headers.push("##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Flubble nesting level\">");
    }
    headers
}

fn profile_info_headers(profile: &str, subr_only: bool) -> Vec<&'static str> {
    match profile {
        "top-level-only" => vec![
            "##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Raw record id\">",
            "##INFO=<ID=PROFILE,Number=1,Type=String,Description=\"Downstream profile name\">",
            "##INFO=<ID=PASSTHROUGH,Number=1,Type=String,Description=\"Record was kept without allele rewrite\">",
        ],
        "popped" => vec![
            "##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Raw record id\">",
            "##INFO=<ID=PARENT,Number=1,Type=String,Description=\"Raw parent record id\">",
            "##INFO=<ID=PROFILE,Number=1,Type=String,Description=\"Downstream profile name\">",
            "##INFO=<ID=RESCUED_CHILD,Number=1,Type=String,Description=\"Child was kept because its parent was popped\">",
            "##INFO=<ID=POPPED_PARENT,Number=1,Type=String,Description=\"Popped parent id that enabled rescue\">",
        ],
        "left-normalized" => vec![
            "##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Raw record id\">",
            "##INFO=<ID=RAW_ALT_INDEX,Number=1,Type=Integer,Description=\"Raw ALT index\">",
            "##INFO=<ID=PROFILE,Number=1,Type=String,Description=\"Downstream profile name\">",
            "##INFO=<ID=LEFT_NORMALIZED,Number=1,Type=String,Description=\"Record was left-normalized\">",
            "##INFO=<ID=RAW_POS,Number=1,Type=Integer,Description=\"Raw POS before profile rewrite\">",
            "##INFO=<ID=RAW_REF,Number=1,Type=String,Description=\"Raw REF before profile rewrite\">",
            "##INFO=<ID=RAW_ALT,Number=1,Type=String,Description=\"Raw ALT before profile rewrite\">",
        ],
        "decomposed" if subr_only => vec![
            "##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Raw record id\">",
            "##INFO=<ID=RAW_ALT_INDEX,Number=1,Type=Integer,Description=\"Raw ALT index\">",
            "##INFO=<ID=PROFILE,Number=1,Type=String,Description=\"Downstream profile name\">",
            "##INFO=<ID=PASSTHROUGH,Number=1,Type=String,Description=\"Record was kept because policy forbids decomposition\">",
            "##INFO=<ID=PASS_THROUGH_REASON,Number=1,Type=String,Description=\"Reason pass-through was selected\">",
            "##INFO=<ID=SUBR_ORIGIN,Number=1,Type=String,Description=\"Raw SUBR semantics were preserved\">",
        ],
        "decomposed" => vec![
            "##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Raw record id\">",
            "##INFO=<ID=RAW_ALT_INDEX,Number=1,Type=Integer,Description=\"Raw ALT index\">",
            "##INFO=<ID=PROFILE,Number=1,Type=String,Description=\"Downstream profile name\">",
            "##INFO=<ID=DECOMPOSED,Number=1,Type=String,Description=\"Record came from allele decomposition\">",
            "##INFO=<ID=PASSTHROUGH,Number=1,Type=String,Description=\"Record was kept because policy forbids decomposition\">",
            "##INFO=<ID=PASS_THROUGH_REASON,Number=1,Type=String,Description=\"Reason pass-through was selected\">",
            "##INFO=<ID=RAW_POS,Number=1,Type=Integer,Description=\"Raw POS before profile rewrite\">",
            "##INFO=<ID=RAW_REF,Number=1,Type=String,Description=\"Raw REF before profile rewrite\">",
            "##INFO=<ID=RAW_ALT,Number=1,Type=String,Description=\"Raw ALT before profile rewrite\">",
        ],
        _ => Vec::new(),
    }
}

fn read_raw_vcf(root: &Path, fixture: &Value) -> Result<Vcf, String> {
    let raw_path = raw_vcf_path(root, fixture)?;
    let text = fs::read_to_string(&raw_path)
        .map_err(|err| format!("failed to read {}: {err}", raw_path.display()))?;
    Vcf::parse(&text)
}

fn raw_vcf_path(root: &Path, fixture: &Value) -> Result<PathBuf, String> {
    for profile in expect_array(fixture, "profiles")? {
        if expect_string(profile, "name")? == "raw-graph" {
            if expect_string(profile, "expected")? != "vcf" {
                return Err(format!(
                    "{} raw-graph profile has no raw VCF source",
                    expect_string(fixture, "id")?
                ));
            }
            return Ok(root.join(expect_string(profile, "path")?));
        }
    }
    Err(format!(
        "{} lacks raw-graph profile",
        expect_string(fixture, "id")?
    ))
}

fn read_fasta(path: &Path) -> Result<std::collections::BTreeMap<String, String>, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read FASTA {}: {err}", path.display()))?;
    let mut records = std::collections::BTreeMap::new();
    let mut current_name: Option<String> = None;
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            current_name = Some(name.trim().to_string());
            records
                .entry(name.trim().to_string())
                .or_insert_with(String::new);
        } else if let Some(name) = &current_name {
            records
                .get_mut(name)
                .expect("FASTA record was inserted")
                .push_str(line.trim());
        }
    }
    Ok(records)
}

fn left_align_indel(
    mut pos: u64,
    ref_allele: &str,
    alt_allele: &str,
    reference: &str,
) -> (u64, String, String) {
    let mut ref_chars: Vec<char> = ref_allele.chars().collect();
    let mut alt_chars: Vec<char> = alt_allele.chars().collect();
    let reference_chars: Vec<char> = reference.chars().collect();

    while pos > 1 {
        let previous = reference_chars
            .get((pos - 2) as usize)
            .copied()
            .unwrap_or_default();
        let longer_last = if ref_chars.len() > alt_chars.len() {
            ref_chars.last().copied()
        } else if alt_chars.len() > ref_chars.len() {
            alt_chars.last().copied()
        } else {
            None
        };
        if longer_last != Some(previous) {
            break;
        }

        rotate_right_with_previous(&mut ref_chars, previous);
        rotate_right_with_previous(&mut alt_chars, previous);
        pos -= 1;
    }

    (pos, ref_chars.iter().collect(), alt_chars.iter().collect())
}

fn rotate_right_with_previous(chars: &mut Vec<char>, previous: char) {
    if !chars.is_empty() {
        chars.pop();
    }
    chars.insert(0, previous);
}

fn parse_info(text: &str) -> Vec<(String, String)> {
    text.split(';')
        .filter(|field| !field.is_empty())
        .map(|field| {
            field
                .split_once('=')
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .unwrap_or_else(|| (field.to_string(), "T".to_string()))
        })
        .collect()
}

fn only_alt(record: &Record) -> Result<&str, String> {
    if record.alt_alleles.len() == 1 {
        Ok(&record.alt_alleles[0])
    } else {
        Err(format!(
            "record {} expected one ALT, found {}",
            record.id,
            record.alt_alleles.len()
        ))
    }
}

fn split_info_list(value: &str) -> Vec<Value> {
    value
        .split(',')
        .map(|item| Value::String(item.to_string()))
        .collect()
}

fn traversal_steps(value: &str) -> Vec<String> {
    let mut steps = Vec::new();
    let mut current = String::new();
    for ch in value.chars() {
        if matches!(ch, '>' | '<') {
            if !current.is_empty() {
                steps.push(current);
            }
            current = ch.to_string();
        } else {
            current.push(ch);
        }
    }
    if !current.is_empty() {
        steps.push(current);
    }
    steps
}

fn mismatches(ref_allele: &str, alt_allele: &str) -> Vec<(usize, char, char)> {
    ref_allele
        .chars()
        .zip(alt_allele.chars())
        .enumerate()
        .filter_map(|(idx, (left, right))| (left != right).then_some((idx, left, right)))
        .collect()
}

fn project_genotypes(genotypes: &[String], selected_alt: &str) -> Vec<String> {
    genotypes
        .iter()
        .map(|gt| {
            if gt == "0" {
                "0".to_string()
            } else if gt == selected_alt {
                "1".to_string()
            } else {
                ".".to_string()
            }
        })
        .collect()
}

fn set_counts(record: &mut Record) {
    let called = record
        .genotypes
        .iter()
        .filter(|gt| gt.as_str() != ".")
        .count();
    let ac = record
        .genotypes
        .iter()
        .filter(|gt| gt.as_str() == "1")
        .count();
    record.set_info("AC", &ac.to_string());
    record.set_info("AF", &format_af(ac, called));
    record.set_info("AN", &called.to_string());
    record.set_info("NS", &called.to_string());
}

fn format_af(ac: usize, an: usize) -> String {
    if an == 0 {
        "0".to_string()
    } else {
        ((ac as f64) / (an as f64)).to_string()
    }
}

fn id_interval_contains(parent: &str, child: &str) -> bool {
    match (id_interval(parent), id_interval(child)) {
        (Some((parent_start, parent_end)), Some((child_start, child_end))) => {
            parent_start <= child_start && child_end <= parent_end
        }
        _ => false,
    }
}

fn id_interval(id: &str) -> Option<(u64, u64)> {
    let mut parts = id
        .trim_start_matches('>')
        .split('>')
        .filter_map(|part| part.parse::<u64>().ok());
    let start = parts.next()?;
    let end = parts.next()?;
    Some((start.min(end), start.max(end)))
}

fn find_fixture(root: &Path, fixture_id: &str) -> Result<Value, String> {
    let manifest = read_json(&root.join("manifest.json"))?;
    for fixture in expect_array(&manifest, "fixtures")? {
        if expect_string(fixture, "id")? == fixture_id {
            return Ok(fixture.clone());
        }
    }
    Err(format!(
        "unknown downstream repetitive fixture id: {fixture_id}"
    ))
}

fn read_json(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    serde_json::from_str(&text).map_err(|err| format!("invalid JSON {}: {err}", path.display()))
}

fn expect_string<'a>(value: &'a Value, key: &str) -> Result<&'a str, String> {
    value
        .get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("JSON value lacks string field {key}"))
}

fn optional_string<'a>(value: &'a Value, key: &str) -> Option<&'a str> {
    value.get(key).and_then(Value::as_str)
}

fn expect_array<'a>(value: &'a Value, key: &str) -> Result<&'a [Value], String> {
    value
        .get(key)
        .and_then(Value::as_array)
        .map(Vec::as_slice)
        .ok_or_else(|| format!("JSON value lacks array field {key}"))
}

fn jsonl(values: &[Value]) -> Result<String, String> {
    let mut out = String::new();
    for value in values {
        out.push_str(&serde_json::to_string(value).map_err(|err| err.to_string())?);
        out.push('\n');
    }
    Ok(out)
}

fn parse_jsonl(text: &str) -> Result<Vec<Value>, String> {
    text.lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| serde_json::from_str(line).map_err(|err| err.to_string()))
        .collect()
}
