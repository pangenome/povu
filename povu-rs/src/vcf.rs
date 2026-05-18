//! Semantic VCF model and deterministic writer for the Rust API.
//!
//! This module is intentionally positioned at the Lean `VariantCall` boundary:
//! callers must provide already verified graph variant sources, allele
//! spellings, traversal strings, allele counts, and genotype columns. The writer
//! validates those semantic obligations before producing VCF text.

use crate::{Error, Result};
use std::collections::HashSet;
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfDocument {
    pub source: String,
    pub samples: Vec<String>,
    pub contigs: Vec<Contig>,
    pub calls: Vec<VariantCall>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub id: String,
    pub length: Option<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantSource {
    Flubble {
        site_id: String,
        parent_id: Option<String>,
        level: u32,
        order_primary: u64,
        order_secondary: u64,
    },
    Hairpin {
        boundary_id: String,
        order_primary: u64,
        order_secondary: u64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    Del,
    Ins,
    Sub,
    Subr,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlleleConstruction {
    Deletion {
        anchor: String,
        deleted: String,
    },
    Insertion {
        anchor: String,
        inserted: String,
    },
    Substitution {
        reference: String,
        alternate: String,
    },
    ReverseSubstitution {
        reference: String,
        alternate: String,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlternateAllele {
    pub construction: AlleleConstruction,
    pub traversal: String,
    pub count: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GenotypeAllele {
    Missing,
    Ref,
    Alt(usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenotypeColumn {
    pub sample: String,
    pub alleles: Vec<GenotypeAllele>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantCall {
    pub source: VariantSource,
    pub chrom: String,
    pub contig_order: u32,
    pub pos: u64,
    pub id: String,
    pub ref_allele: String,
    pub ref_traversal: String,
    pub alternates: Vec<AlternateAllele>,
    pub variant_type: VariantType,
    pub tangled: bool,
    pub reference_allele_count: u64,
    pub genotypes: Vec<GenotypeColumn>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Info {
    pub ac: Vec<u64>,
    pub af: Vec<AlleleFrequency>,
    pub an: u64,
    pub ns: usize,
    pub at: Vec<String>,
    pub variant_type: VariantType,
    pub tangled: bool,
    pub enclosing_site: Option<String>,
    pub level: Option<u32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlleleFrequency {
    pub alternate_count: u64,
    pub total_alleles: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub source: VariantSource,
    pub chrom: String,
    pub contig_order: u32,
    pub pos: u64,
    pub id: String,
    pub ref_allele: String,
    pub alternates: Vec<AlternateAllele>,
    pub qual: String,
    pub filter: String,
    pub info: Info,
    pub format: String,
    pub reference_allele_count: u64,
    pub ref_traversal: String,
    pub genotypes: Vec<GenotypeColumn>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct OrderKey {
    pub contig_order: u32,
    pub pos: u64,
    pub source_primary: u64,
    pub source_secondary: u64,
}

impl VcfDocument {
    pub fn new<S, I>(samples: I, calls: Vec<VariantCall>) -> Self
    where
        S: Into<String>,
        I: IntoIterator<Item = S>,
    {
        Self {
            source: "povu-rs".to_string(),
            samples: samples.into_iter().map(Into::into).collect(),
            contigs: Vec::new(),
            calls,
        }
    }

    pub fn with_source(mut self, source: impl Into<String>) -> Self {
        self.source = source.into();
        self
    }

    pub fn with_contigs(mut self, contigs: Vec<Contig>) -> Self {
        self.contigs = contigs;
        self
    }

    pub fn records(&self) -> Result<Vec<Record>> {
        self.validate_metadata()?;

        let mut ids = HashSet::new();
        let mut records = Vec::with_capacity(self.calls.len());
        for call in &self.calls {
            self.validate_call_samples(call)?;

            if !ids.insert(call.id.as_str()) {
                return Err(Error::invalid_vcf(format!(
                    "duplicate VCF record id '{}'",
                    call.id
                )));
            }

            records.push(Record::from_call(call.clone())?);
        }

        records.sort_by_key(Record::order_key);
        Ok(records)
    }

    pub fn to_vcf_string(&self) -> Result<String> {
        let records = self.records()?;
        let mut out = String::new();

        push_line(&mut out, "##fileformat=VCFv4.2");
        push_line(&mut out, &format!("##source={}", self.source));
        for contig in &self.contigs {
            match contig.length {
                Some(length) => push_line(
                    &mut out,
                    &format!("##contig=<ID={},length={}>", contig.id, length),
                ),
                None => push_line(&mut out, &format!("##contig=<ID={}>", contig.id)),
            }
        }
        push_line(
            &mut out,
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate allele count\">",
        );
        push_line(
            &mut out,
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency\">",
        );
        push_line(
            &mut out,
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total called allele count\">",
        );
        push_line(
            &mut out,
            "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">",
        );
        push_line(
            &mut out,
            "##INFO=<ID=AT,Number=R,Type=String,Description=\"Graph traversal provenance for REF and ALT alleles\">",
        );
        push_line(
            &mut out,
            "##INFO=<ID=VARTYPE,Number=1,Type=String,Description=\"povu variant type\">",
        );
        push_line(
            &mut out,
            "##INFO=<ID=TANGLED,Number=1,Type=String,Description=\"Graph-faithful complex representation marker\">",
        );
        push_line(
            &mut out,
            "##INFO=<ID=ES,Number=1,Type=String,Description=\"Enclosing flubble site id\">",
        );
        push_line(
            &mut out,
            "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Flubble nesting level\">",
        );
        push_line(
            &mut out,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        );

        let mut header = vec![
            "#CHROM".to_string(),
            "POS".to_string(),
            "ID".to_string(),
            "REF".to_string(),
            "ALT".to_string(),
            "QUAL".to_string(),
            "FILTER".to_string(),
            "INFO".to_string(),
            "FORMAT".to_string(),
        ];
        header.extend(self.samples.iter().cloned());
        push_line(&mut out, &join_tab(header));

        for record in &records {
            push_line(&mut out, &record.to_vcf_line());
        }

        Ok(out)
    }

    pub fn write_path(&self, path: impl AsRef<Path>) -> Result<()> {
        fs::write(path, self.to_vcf_string()?)?;
        Ok(())
    }

    fn validate_metadata(&self) -> Result<()> {
        if self.source.is_empty() {
            return Err(Error::invalid_vcf("VCF source must be nonempty"));
        }
        if self.samples.is_empty() {
            return Err(Error::invalid_vcf(
                "VCF output requires at least one sample column",
            ));
        }

        let mut samples = HashSet::new();
        for sample in &self.samples {
            if sample.is_empty() {
                return Err(Error::invalid_vcf("sample names must be nonempty"));
            }
            if !samples.insert(sample.as_str()) {
                return Err(Error::invalid_vcf(format!(
                    "duplicate sample name '{}'",
                    sample
                )));
            }
        }

        let mut contigs = HashSet::new();
        for contig in &self.contigs {
            if contig.id.is_empty() {
                return Err(Error::invalid_vcf("contig ids must be nonempty"));
            }
            if !contigs.insert(contig.id.as_str()) {
                return Err(Error::invalid_vcf(format!(
                    "duplicate contig id '{}'",
                    contig.id
                )));
            }
        }

        Ok(())
    }

    fn validate_call_samples(&self, call: &VariantCall) -> Result<()> {
        if call.genotypes.len() != self.samples.len() {
            return Err(Error::invalid_vcf(format!(
                "record '{}' has {} genotype columns but document has {} samples",
                call.id,
                call.genotypes.len(),
                self.samples.len()
            )));
        }

        for (expected, column) in self.samples.iter().zip(&call.genotypes) {
            if &column.sample != expected {
                return Err(Error::invalid_vcf(format!(
                    "record '{}' genotype sample '{}' does not match header sample '{}'",
                    call.id, column.sample, expected
                )));
            }
        }

        Ok(())
    }
}

impl Contig {
    pub fn new(id: impl Into<String>, length: Option<u64>) -> Self {
        Self {
            id: id.into(),
            length,
        }
    }
}

impl VariantSource {
    pub fn flubble(
        site_id: impl Into<String>,
        level: u32,
        order_primary: u64,
        order_secondary: u64,
    ) -> Self {
        Self::Flubble {
            site_id: site_id.into(),
            parent_id: None,
            level,
            order_primary,
            order_secondary,
        }
    }

    pub fn nested_flubble(
        site_id: impl Into<String>,
        parent_id: impl Into<String>,
        level: u32,
        order_primary: u64,
        order_secondary: u64,
    ) -> Self {
        Self::Flubble {
            site_id: site_id.into(),
            parent_id: Some(parent_id.into()),
            level,
            order_primary,
            order_secondary,
        }
    }

    pub fn hairpin(
        boundary_id: impl Into<String>,
        order_primary: u64,
        order_secondary: u64,
    ) -> Self {
        Self::Hairpin {
            boundary_id: boundary_id.into(),
            order_primary,
            order_secondary,
        }
    }

    pub fn order_primary(&self) -> u64 {
        match self {
            VariantSource::Flubble { order_primary, .. }
            | VariantSource::Hairpin { order_primary, .. } => *order_primary,
        }
    }

    pub fn order_secondary(&self) -> u64 {
        match self {
            VariantSource::Flubble {
                order_secondary, ..
            }
            | VariantSource::Hairpin {
                order_secondary, ..
            } => *order_secondary,
        }
    }

    pub fn supports_variant_type(&self, variant_type: VariantType) -> bool {
        match self {
            VariantSource::Flubble { .. } => variant_type.is_flubble_type(),
            VariantSource::Hairpin { .. } => variant_type == VariantType::Subr,
        }
    }

    fn enclosing_site(&self) -> Option<&str> {
        match self {
            VariantSource::Flubble { site_id, .. } => Some(site_id.as_str()),
            VariantSource::Hairpin { .. } => None,
        }
    }

    fn level(&self) -> Option<u32> {
        match self {
            VariantSource::Flubble { level, .. } => Some(*level),
            VariantSource::Hairpin { .. } => None,
        }
    }

    fn validate(&self) -> Result<()> {
        match self {
            VariantSource::Flubble { site_id, .. } if site_id.is_empty() => Err(
                Error::invalid_vcf("flubble source site id must be nonempty"),
            ),
            VariantSource::Hairpin { boundary_id, .. } if boundary_id.is_empty() => Err(
                Error::invalid_vcf("hairpin source boundary id must be nonempty"),
            ),
            _ => Ok(()),
        }
    }
}

impl VariantType {
    pub fn vcf_label(self) -> &'static str {
        match self {
            VariantType::Del => "DEL",
            VariantType::Ins => "INS",
            VariantType::Sub => "SUB",
            VariantType::Subr => "SUBR",
        }
    }

    pub fn is_flubble_type(self) -> bool {
        matches!(self, VariantType::Del | VariantType::Ins | VariantType::Sub)
    }
}

impl AlleleConstruction {
    pub fn deletion(anchor: impl Into<String>, deleted: impl Into<String>) -> Self {
        Self::Deletion {
            anchor: anchor.into(),
            deleted: deleted.into(),
        }
    }

    pub fn insertion(anchor: impl Into<String>, inserted: impl Into<String>) -> Self {
        Self::Insertion {
            anchor: anchor.into(),
            inserted: inserted.into(),
        }
    }

    pub fn substitution(reference: impl Into<String>, alternate: impl Into<String>) -> Self {
        Self::Substitution {
            reference: reference.into(),
            alternate: alternate.into(),
        }
    }

    pub fn reverse_substitution(
        reference: impl Into<String>,
        alternate: impl Into<String>,
    ) -> Self {
        Self::ReverseSubstitution {
            reference: reference.into(),
            alternate: alternate.into(),
        }
    }

    pub fn reference_allele(&self) -> String {
        match self {
            AlleleConstruction::Deletion { anchor, deleted } => {
                format!("{}{}", anchor, deleted)
            }
            AlleleConstruction::Insertion { anchor, .. } => anchor.clone(),
            AlleleConstruction::Substitution { reference, .. }
            | AlleleConstruction::ReverseSubstitution { reference, .. } => reference.clone(),
        }
    }

    pub fn alternate_allele(&self) -> String {
        match self {
            AlleleConstruction::Deletion { anchor, .. } => anchor.clone(),
            AlleleConstruction::Insertion { anchor, inserted } => {
                format!("{}{}", anchor, inserted)
            }
            AlleleConstruction::Substitution { alternate, .. }
            | AlleleConstruction::ReverseSubstitution { alternate, .. } => alternate.clone(),
        }
    }

    pub fn variant_type(&self) -> VariantType {
        match self {
            AlleleConstruction::Deletion { .. } => VariantType::Del,
            AlleleConstruction::Insertion { .. } => VariantType::Ins,
            AlleleConstruction::Substitution { .. } => VariantType::Sub,
            AlleleConstruction::ReverseSubstitution { .. } => VariantType::Subr,
        }
    }

    fn validate(&self) -> Result<()> {
        match self {
            AlleleConstruction::Deletion { anchor, deleted } => {
                if anchor.is_empty() || deleted.is_empty() {
                    return Err(Error::invalid_vcf(
                        "deletion construction requires nonempty anchor and deleted sequence",
                    ));
                }
            }
            AlleleConstruction::Insertion { anchor, inserted } => {
                if anchor.is_empty() || inserted.is_empty() {
                    return Err(Error::invalid_vcf(
                        "insertion construction requires nonempty anchor and inserted sequence",
                    ));
                }
            }
            AlleleConstruction::Substitution {
                reference,
                alternate,
            }
            | AlleleConstruction::ReverseSubstitution {
                reference,
                alternate,
            } => {
                if reference.is_empty() || alternate.is_empty() {
                    return Err(Error::invalid_vcf(
                        "substitution construction requires nonempty reference and alternate",
                    ));
                }
            }
        }

        if self.reference_allele().is_empty() || self.alternate_allele().is_empty() {
            return Err(Error::invalid_vcf(
                "constructed REF and ALT alleles must be nonempty",
            ));
        }

        Ok(())
    }
}

impl AlternateAllele {
    pub fn new(construction: AlleleConstruction, traversal: impl Into<String>, count: u64) -> Self {
        Self {
            construction,
            traversal: traversal.into(),
            count,
        }
    }

    pub fn reference_allele(&self) -> String {
        self.construction.reference_allele()
    }

    pub fn alternate_allele(&self) -> String {
        self.construction.alternate_allele()
    }

    pub fn variant_type(&self) -> VariantType {
        self.construction.variant_type()
    }

    fn validate(&self) -> Result<()> {
        self.construction.validate()?;
        if self.traversal.is_empty() {
            return Err(Error::invalid_vcf("ALT traversal strings must be nonempty"));
        }
        Ok(())
    }
}

impl GenotypeColumn {
    pub fn new(sample: impl Into<String>, alleles: Vec<GenotypeAllele>) -> Self {
        Self {
            sample: sample.into(),
            alleles,
        }
    }

    pub fn has_data(&self) -> bool {
        self.alleles.iter().any(GenotypeAllele::is_called)
    }

    fn validate(&self, alt_count: usize) -> Result<()> {
        if self.sample.is_empty() {
            return Err(Error::invalid_vcf("genotype sample names must be nonempty"));
        }
        if self.alleles.is_empty() {
            return Err(Error::invalid_vcf(format!(
                "genotype column '{}' must contain at least one phase",
                self.sample
            )));
        }

        for allele in &self.alleles {
            if let GenotypeAllele::Alt(index) = allele {
                if *index == 0 || *index > alt_count {
                    return Err(Error::invalid_vcf(format!(
                        "genotype column '{}' references missing ALT index {}",
                        self.sample, index
                    )));
                }
            }
        }

        Ok(())
    }
}

impl GenotypeAllele {
    pub fn is_called(&self) -> bool {
        !matches!(self, GenotypeAllele::Missing)
    }

    fn to_vcf_token(&self) -> String {
        match self {
            GenotypeAllele::Missing => ".".to_string(),
            GenotypeAllele::Ref => "0".to_string(),
            GenotypeAllele::Alt(index) => index.to_string(),
        }
    }
}

impl VariantCall {
    pub fn validate(&self) -> Result<()> {
        if self.chrom.is_empty() {
            return Err(Error::invalid_vcf("CHROM must be nonempty"));
        }
        if self.pos == 0 {
            return Err(Error::invalid_vcf("POS must be one-based and positive"));
        }
        if self.id.is_empty() {
            return Err(Error::invalid_vcf("ID must be nonempty"));
        }
        if self.ref_allele.is_empty() {
            return Err(Error::invalid_vcf("REF must be nonempty"));
        }
        if self.ref_traversal.is_empty() {
            return Err(Error::invalid_vcf("REF traversal string must be nonempty"));
        }
        if self.alternates.is_empty() {
            return Err(Error::invalid_vcf(
                "records require at least one ALT allele",
            ));
        }

        self.source.validate()?;
        if !self.source.supports_variant_type(self.variant_type) {
            return Err(Error::invalid_vcf(format!(
                "source does not support variant type {}",
                self.variant_type.vcf_label()
            )));
        }

        for alt in &self.alternates {
            alt.validate()?;
            if alt.reference_allele() != self.ref_allele {
                return Err(Error::invalid_vcf(format!(
                    "ALT construction REF '{}' does not match record REF '{}'",
                    alt.reference_allele(),
                    self.ref_allele
                )));
            }
            if alt.variant_type() != self.variant_type {
                return Err(Error::invalid_vcf(format!(
                    "ALT construction type {} does not match record VARTYPE {}",
                    alt.variant_type().vcf_label(),
                    self.variant_type.vcf_label()
                )));
            }
        }

        if self.total_alleles() == 0 {
            return Err(Error::invalid_vcf("total allele count must be positive"));
        }

        for column in &self.genotypes {
            column.validate(self.alternates.len())?;
        }

        Ok(())
    }

    pub fn alternate_counts(&self) -> Vec<u64> {
        self.alternates.iter().map(|alt| alt.count).collect()
    }

    pub fn total_alleles(&self) -> u64 {
        self.reference_allele_count + self.alternates.iter().map(|alt| alt.count).sum::<u64>()
    }

    pub fn order_key(&self) -> OrderKey {
        OrderKey {
            contig_order: self.contig_order,
            pos: self.pos,
            source_primary: self.source.order_primary(),
            source_secondary: self.source.order_secondary(),
        }
    }
}

impl Record {
    pub fn from_call(call: VariantCall) -> Result<Self> {
        call.validate()?;

        let an = call.total_alleles();
        let ac = call.alternate_counts();
        let af = ac
            .iter()
            .map(|count| AlleleFrequency {
                alternate_count: *count,
                total_alleles: an,
            })
            .collect();
        let mut at = Vec::with_capacity(1 + call.alternates.len());
        at.push(call.ref_traversal.clone());
        at.extend(call.alternates.iter().map(|alt| alt.traversal.clone()));

        let info = Info {
            ac,
            af,
            an,
            ns: call
                .genotypes
                .iter()
                .filter(|column| column.has_data())
                .count(),
            at,
            variant_type: call.variant_type,
            tangled: call.tangled,
            enclosing_site: call.source.enclosing_site().map(str::to_string),
            level: call.source.level(),
        };

        Ok(Self {
            source: call.source,
            chrom: call.chrom,
            contig_order: call.contig_order,
            pos: call.pos,
            id: call.id,
            ref_allele: call.ref_allele,
            alternates: call.alternates,
            qual: "60".to_string(),
            filter: "PASS".to_string(),
            info,
            format: "GT".to_string(),
            reference_allele_count: call.reference_allele_count,
            ref_traversal: call.ref_traversal,
            genotypes: call.genotypes,
        })
    }

    pub fn alternate_alleles(&self) -> Vec<String> {
        self.alternates
            .iter()
            .map(AlternateAllele::alternate_allele)
            .collect()
    }

    pub fn order_key(&self) -> OrderKey {
        OrderKey {
            contig_order: self.contig_order,
            pos: self.pos,
            source_primary: self.source.order_primary(),
            source_secondary: self.source.order_secondary(),
        }
    }

    fn to_vcf_line(&self) -> String {
        let mut fields = vec![
            self.chrom.clone(),
            self.pos.to_string(),
            self.id.clone(),
            self.ref_allele.clone(),
            join_comma(self.alternate_alleles()),
            self.qual.clone(),
            self.filter.clone(),
            self.info.to_vcf_string(),
            self.format.clone(),
        ];
        fields.extend(self.genotypes.iter().map(format_genotype_column));
        join_tab(fields)
    }
}

impl Info {
    fn to_vcf_string(&self) -> String {
        let mut fields = vec![
            format!("AC={}", join_numbers(&self.ac)),
            format!(
                "AF={}",
                join_comma(self.af.iter().map(format_allele_frequency))
            ),
            format!("AN={}", self.an),
            format!("NS={}", self.ns),
            format!("AT={}", join_comma(self.at.iter().cloned())),
            format!("VARTYPE={}", self.variant_type.vcf_label()),
            format!("TANGLED={}", if self.tangled { "T" } else { "F" }),
        ];

        if let Some(site) = &self.enclosing_site {
            fields.push(format!("ES={}", site));
        }
        if let Some(level) = self.level {
            fields.push(format!("LV={}", level));
        }

        fields.join(";")
    }
}

fn format_genotype_column(column: &GenotypeColumn) -> String {
    if column
        .alleles
        .iter()
        .all(|allele| matches!(allele, GenotypeAllele::Missing))
    {
        return ".".to_string();
    }

    column
        .alleles
        .iter()
        .map(GenotypeAllele::to_vcf_token)
        .collect::<Vec<_>>()
        .join("|")
}

fn format_allele_frequency(freq: &AlleleFrequency) -> String {
    if freq.total_alleles == 0 {
        return "NaN".to_string();
    }
    ((freq.alternate_count as f64) / (freq.total_alleles as f64)).to_string()
}

fn join_numbers(values: &[u64]) -> String {
    values
        .iter()
        .map(u64::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn join_comma<I, S>(values: I) -> String
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    values
        .into_iter()
        .map(Into::into)
        .collect::<Vec<_>>()
        .join(",")
}

fn join_tab<I, S>(values: I) -> String
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    values
        .into_iter()
        .map(Into::into)
        .collect::<Vec<_>>()
        .join("\t")
}

fn push_line(out: &mut String, line: &str) {
    out.push_str(line);
    out.push('\n');
}
