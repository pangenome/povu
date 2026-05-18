use povu::{
    AlleleConstruction, AlternateAllele, Contig, Error, GenotypeAllele, GenotypeColumn,
    VariantCall, VariantSource, VariantType, VcfDocument,
};

fn ref_sample(sample: &str) -> GenotypeColumn {
    GenotypeColumn::new(sample, vec![GenotypeAllele::Ref])
}

fn alt_sample(sample: &str) -> GenotypeColumn {
    GenotypeColumn::new(sample, vec![GenotypeAllele::Alt(1)])
}

fn missing_sample(sample: &str) -> GenotypeColumn {
    GenotypeColumn::new(sample, vec![GenotypeAllele::Missing])
}

fn flubble_call(
    id: &str,
    pos: u64,
    ref_allele: &str,
    ref_traversal: &str,
    alt: AlternateAllele,
    variant_type: VariantType,
    level: u32,
    order_primary: u64,
    order_secondary: u64,
    genotypes: Vec<GenotypeColumn>,
) -> VariantCall {
    VariantCall {
        source: VariantSource::flubble(id, level, order_primary, order_secondary),
        chrom: "HG1#1#chr1".to_string(),
        contig_order: 0,
        pos,
        id: id.to_string(),
        ref_allele: ref_allele.to_string(),
        ref_traversal: ref_traversal.to_string(),
        alternates: vec![alt],
        variant_type,
        tangled: false,
        reference_allele_count: 1,
        genotypes,
    }
}

fn minimal_substitution_call() -> VariantCall {
    flubble_call(
        ">0>3",
        2,
        "C",
        ">1",
        AlternateAllele::new(AlleleConstruction::substitution("C", "G"), ">2", 1),
        VariantType::Sub,
        0,
        0,
        3,
        vec![ref_sample("HG1"), alt_sample("HG2")],
    )
}

fn insertion_flubble_call() -> VariantCall {
    flubble_call(
        ">0>1",
        1,
        "A",
        ">0",
        AlternateAllele::new(AlleleConstruction::insertion("A", "G"), ">0>2", 1),
        VariantType::Ins,
        0,
        0,
        1,
        vec![ref_sample("HG1"), alt_sample("HG2")],
    )
}

fn deletion_flubble_call() -> VariantCall {
    flubble_call(
        ">0>1",
        1,
        "AG",
        ">0>2",
        AlternateAllele::new(AlleleConstruction::deletion("A", "G"), ">0", 1),
        VariantType::Del,
        0,
        0,
        1,
        vec![ref_sample("HG1"), alt_sample("HG2")],
    )
}

fn nested_deletion_call() -> VariantCall {
    VariantCall {
        source: VariantSource::nested_flubble(">1>4", ">0>5", 1, 1, 4),
        chrom: "HG1#1#chr1".to_string(),
        contig_order: 0,
        pos: 2,
        id: ">1>4".to_string(),
        ref_allele: "CT".to_string(),
        ref_traversal: ">1>3".to_string(),
        alternates: vec![AlternateAllele::new(
            AlleleConstruction::deletion("C", "T"),
            ">1",
            1,
        )],
        variant_type: VariantType::Del,
        tangled: false,
        reference_allele_count: 1,
        genotypes: vec![ref_sample("HG1"), alt_sample("HG2"), missing_sample("HG3")],
    }
}

fn hairpin_inversion_subr_call() -> VariantCall {
    VariantCall {
        source: VariantSource::hairpin(">1>5", 1, 5),
        chrom: "ref".to_string(),
        contig_order: 0,
        pos: 2,
        id: ">1>5".to_string(),
        ref_allele: "ACGTA".to_string(),
        ref_traversal: ">1>2>3>4>5".to_string(),
        alternates: vec![AlternateAllele::new(
            AlleleConstruction::reverse_substitution("ACGTA", "TACGT"),
            "<5<4<3<2<1",
            1,
        )],
        variant_type: VariantType::Subr,
        tangled: false,
        reference_allele_count: 1,
        genotypes: vec![ref_sample("ref"), alt_sample("alt")],
    }
}

fn second_ordered_substitution_call() -> VariantCall {
    flubble_call(
        ">3>6",
        4,
        "A",
        ">4",
        AlternateAllele::new(AlleleConstruction::substitution("A", "C"), ">5", 1),
        VariantType::Sub,
        0,
        3,
        6,
        vec![ref_sample("HG1"), alt_sample("HG2")],
    )
}

fn record_lines(text: &str) -> Vec<&str> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .collect::<Vec<_>>()
}

#[test]
fn header_matches_modernized_rust_contract() {
    let doc = VcfDocument::new(["HG1", "HG2"], vec![minimal_substitution_call()])
        .with_contigs(vec![Contig::new("HG1#1#chr1", None)]);
    let text = doc.to_vcf_string().unwrap();

    assert!(text.starts_with("##fileformat=VCFv4.2\n##source=povu-rs\n"));
    assert!(text.contains("##contig=<ID=HG1#1#chr1>\n"));
    assert!(text.contains("##INFO=<ID=AN,Number=1,Type=Integer"));
    assert_eq!(text.matches("##FORMAT=<ID=GT,").count(), 1);
    assert!(!text.contains("##fileDate="));
    assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\n"));
}

#[test]
fn minimal_substitution_fixture_record_matches_spec() {
    let text = VcfDocument::new(["HG1", "HG2"], vec![minimal_substitution_call()])
        .to_vcf_string()
        .unwrap();

    assert_eq!(
        record_lines(&text),
        vec![
            "HG1#1#chr1\t2\t>0>3\tC\tG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1,>2;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0\tGT\t0\t1"
        ]
    );
}

#[test]
fn insertion_flubble_fixture_record_matches_spec() {
    let text = VcfDocument::new(["HG1", "HG2"], vec![insertion_flubble_call()])
        .to_vcf_string()
        .unwrap();

    assert_eq!(
        record_lines(&text),
        vec![
            "HG1#1#chr1\t1\t>0>1\tA\tAG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>0,>0>2;VARTYPE=INS;TANGLED=F;ES=>0>1;LV=0\tGT\t0\t1"
        ]
    );
}

#[test]
fn deletion_flubble_fixture_record_matches_spec() {
    let text = VcfDocument::new(["HG1", "HG2"], vec![deletion_flubble_call()])
        .to_vcf_string()
        .unwrap();

    assert_eq!(
        record_lines(&text),
        vec![
            "HG1#1#chr1\t1\t>0>1\tAG\tA\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>0>2,>0;VARTYPE=DEL;TANGLED=F;ES=>0>1;LV=0\tGT\t0\t1"
        ]
    );
}

#[test]
fn nested_deletion_fixture_record_matches_spec() {
    let text = VcfDocument::new(["HG1", "HG2", "HG3"], vec![nested_deletion_call()])
        .to_vcf_string()
        .unwrap();

    assert_eq!(
        record_lines(&text),
        vec![
            "HG1#1#chr1\t2\t>1>4\tCT\tC\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1>3,>1;VARTYPE=DEL;TANGLED=F;ES=>1>4;LV=1\tGT\t0\t1\t."
        ]
    );
}

#[test]
fn hairpin_inversion_subr_fixture_record_matches_spec() {
    let text = VcfDocument::new(["ref", "alt"], vec![hairpin_inversion_subr_call()])
        .to_vcf_string()
        .unwrap();

    let records = record_lines(&text);
    assert_eq!(
        records,
        vec![
            "ref\t2\t>1>5\tACGTA\tTACGT\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1>2>3>4>5,<5<4<3<2<1;VARTYPE=SUBR;TANGLED=F\tGT\t0\t1"
        ]
    );
    assert!(!records[0].contains(";ES="));
    assert!(!records[0].contains(";LV="));
}

#[test]
fn linear_no_variant_fixture_emits_header_only() {
    let text = VcfDocument::new(["HG1", "HG2"], vec![])
        .to_vcf_string()
        .unwrap();

    assert!(record_lines(&text).is_empty());
    assert!(text.ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\n"));
}

#[test]
fn two_ordered_substitutions_fixture_uses_semantic_order_key() {
    let text = VcfDocument::new(
        ["HG1", "HG2"],
        vec![
            second_ordered_substitution_call(),
            minimal_substitution_call(),
        ],
    )
    .to_vcf_string()
    .unwrap();

    assert_eq!(
        record_lines(&text),
        vec![
            "HG1#1#chr1\t2\t>0>3\tC\tG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1,>2;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0\tGT\t0\t1",
            "HG1#1#chr1\t4\t>3>6\tA\tC\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>4,>5;VARTYPE=SUB;TANGLED=F;ES=>3>6;LV=0\tGT\t0\t1",
        ]
    );
}

#[test]
fn allele_frequency_uses_stable_decimal_precision() {
    let mut call = minimal_substitution_call();
    call.reference_allele_count = 2;
    call.alternates[0].count = 1;
    call.genotypes = vec![ref_sample("S1"), ref_sample("S2"), alt_sample("S3")];

    let text = VcfDocument::new(["S1", "S2", "S3"], vec![call])
        .to_vcf_string()
        .unwrap();

    assert!(record_lines(&text)[0].contains("AF=0.3333333333333333"));
}

#[test]
fn all_missing_phases_collapse_to_single_missing_genotype() {
    let mut call = minimal_substitution_call();
    call.genotypes = vec![GenotypeColumn::new(
        "HG1",
        vec![GenotypeAllele::Missing, GenotypeAllele::Missing],
    )];

    let text = VcfDocument::new(["HG1"], vec![call])
        .to_vcf_string()
        .unwrap();

    assert!(record_lines(&text)[0].ends_with("\tGT\t."));
}

#[test]
fn partial_missing_phases_are_preserved_in_place() {
    let mut call = minimal_substitution_call();
    call.genotypes = vec![GenotypeColumn::new(
        "HG1",
        vec![GenotypeAllele::Ref, GenotypeAllele::Missing],
    )];

    let text = VcfDocument::new(["HG1"], vec![call])
        .to_vcf_string()
        .unwrap();

    assert!(record_lines(&text)[0].ends_with("\tGT\t0|."));
}

#[test]
fn write_path_uses_same_serialization_as_string_output() {
    let doc = VcfDocument::new(["HG1", "HG2"], vec![minimal_substitution_call()]);
    let expected = doc.to_vcf_string().unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.vcf");

    doc.write_path(&path).unwrap();

    assert_eq!(std::fs::read_to_string(path).unwrap(), expected);
}

#[test]
fn invalid_inputs_return_structured_errors() {
    let mut call = hairpin_inversion_subr_call();
    call.source = VariantSource::flubble(">1>5", 0, 1, 5);

    let err = VcfDocument::new(["ref", "alt"], vec![call])
        .to_vcf_string()
        .unwrap_err();

    match err {
        Error::InvalidVcf { message } => {
            assert!(message.contains("source does not support variant type SUBR"));
        }
        other => panic!("expected InvalidVcf, got {other:?}"),
    }
}

#[test]
fn duplicate_record_ids_are_rejected() {
    let err = VcfDocument::new(
        ["HG1", "HG2"],
        vec![minimal_substitution_call(), minimal_substitution_call()],
    )
    .to_vcf_string()
    .unwrap_err();

    match err {
        Error::InvalidVcf { message } => assert!(message.contains("duplicate VCF record id")),
        other => panic!("expected InvalidVcf, got {other:?}"),
    }
}
