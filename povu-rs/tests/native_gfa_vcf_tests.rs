use povu::{
    detect_flubble_stack, gfa_to_vcf, gfa_to_vcf_document, Error, FlubbleCandidate, NativeGfa,
};
use std::fs;
use std::path::{Path, PathBuf};

fn fixture_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("povu-rs has a repository parent")
        .join("tests/lean4_conformance/fixtures")
        .join(name)
}

fn native_record_lines(fixture: &str, reference: &str) -> Vec<String> {
    let temp = tempfile::tempdir().expect("tempdir");
    let ref_file = temp.path().join("refs.txt");
    fs::write(&ref_file, format!("{reference}\n")).expect("write reference file");
    let document = gfa_to_vcf_document(fixture_path(fixture), Some(&ref_file)).expect("native VCF");
    document
        .to_vcf_string()
        .expect("serialize VCF")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

fn native_record_lines_from_text(gfa: &str, reference: &str) -> Vec<String> {
    let temp = tempfile::tempdir().expect("tempdir");
    let gfa_path = temp.path().join("input.gfa");
    let ref_file = temp.path().join("refs.txt");
    fs::write(&gfa_path, gfa).expect("write GFA");
    fs::write(&ref_file, format!("{reference}\n")).expect("write reference file");
    let document = gfa_to_vcf_document(&gfa_path, Some(&ref_file)).expect("native VCF");
    document
        .to_vcf_string()
        .expect("serialize VCF")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

#[test]
fn flubble_stack_port_matches_lean_close_after_gap_rule() {
    let stack = vec![
        FlubbleCandidate { id: 1, class_id: 1 },
        FlubbleCandidate { id: 2, class_id: 2 },
        FlubbleCandidate { id: 3, class_id: 1 },
        FlubbleCandidate { id: 4, class_id: 2 },
    ];
    let boundaries = detect_flubble_stack(&stack);
    assert_eq!(boundaries.len(), 2);
    assert_eq!(boundaries[0].open_edge, 1);
    assert_eq!(boundaries[0].close_edge, 3);
    assert_eq!(boundaries[1].open_edge, 2);
    assert_eq!(boundaries[1].close_edge, 4);

    let immediate_same = vec![
        FlubbleCandidate {
            id: 10,
            class_id: 7,
        },
        FlubbleCandidate {
            id: 11,
            class_id: 7,
        },
        FlubbleCandidate {
            id: 12,
            class_id: 7,
        },
    ];
    assert!(detect_flubble_stack(&immediate_same).is_empty());
}

#[test]
fn native_gfa_to_vcf_matches_lean_minimal_substitution_fixture() {
    assert_eq!(
        native_record_lines("minimal_substitution.gfa", "HG1"),
        vec![
            "HG1#1#chr1\t2\t>0>3\tC\tG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1,>2;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0\tGT\t0\t1"
        ]
    );
}

#[test]
fn native_gfa_to_vcf_parses_w_lines_with_same_path_model() {
    let gfa = "\
H\tVN:Z:1.1
S\t0\tA
S\t1\tC
S\t2\tG
S\t3\tT
L\t0\t+\t1\t+\t0M
L\t1\t+\t3\t+\t0M
L\t0\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
W\tHG1\t1\tchr1\t0\t3\t>0>1>3
W\tHG2\t1\tchr1\t0\t3\t>0>2>3
";
    assert_eq!(
        native_record_lines_from_text(gfa, "HG1"),
        vec![
            "HG1#1#chr1:0-3\t2\t>0>3\tC\tG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1,>2;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0\tGT\t0\t1"
        ]
    );
}

#[test]
fn native_gfa_exposes_flubble_decomposition_sites() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tAA
S\t6\tA
S\t7\tCCC
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t2\t+\t5\t+\t0M
L\t5\t+\t4\t+\t0M
L\t4\t+\t6\t+\t0M
L\t1\t+\t7\t+\t0M
L\t7\t+\t6\t+\t0M
P\tref\t1+,2+,3+,4+,6+\t*
P\tinner_alt\t1+,2+,5+,4+,6+\t*
P\touter_alt\t1+,7+,6+\t*
";
    let graph = NativeGfa::parse(gfa).expect("parse native GFA");
    let decomposition = graph
        .decompose_flubbles(&["ref".to_string()])
        .expect("decompose flubbles");

    assert_eq!(decomposition.reference_path, "ref");
    assert_eq!(
        decomposition
            .sites
            .iter()
            .map(|site| (
                site.id.as_str(),
                site.parent_id.as_deref(),
                site.level,
                site.reference_start_step,
                site.reference_end_step,
                site.is_leaf,
            ))
            .collect::<Vec<_>>(),
        vec![
            (">1>6", None, 0, 0, 4, false),
            (">2>4", Some(">1>6"), 1, 1, 3, true),
        ]
    );

    let leaves = decomposition.leaf_sites_bottom_up();
    assert_eq!(leaves.len(), 1);
    assert_eq!(leaves[0].id, ">2>4");
    assert_eq!(leaves[0].start.token(), ">2");
    assert_eq!(leaves[0].end.token(), ">4");
}

#[test]
fn one_shot_gfa_to_vcf_uses_native_rust_path() {
    let temp = tempfile::tempdir().expect("tempdir");
    let ref_file = temp.path().join("refs.txt");
    let out_vcf = temp.path().join("out.vcf");
    fs::write(&ref_file, "HG1\n").expect("write reference file");

    gfa_to_vcf(
        fixture_path("minimal_substitution.gfa"),
        &out_vcf,
        Some(&ref_file),
    )
    .expect("native one-shot VCF");

    let text = fs::read_to_string(out_vcf).expect("read VCF");
    assert!(text.contains("##source=povu-rs-native\n"));
    assert!(text.contains("HG1#1#chr1\t2\t>0>3\tC\tG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2"));
}

#[test]
fn native_gfa_to_vcf_matches_lean_insertion_and_deletion_fixtures() {
    assert_eq!(
        native_record_lines("insertion_flubble.gfa", "HG1"),
        vec![
            "HG1#1#chr1\t1\t>0>1\tA\tAG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>0,>0>2;VARTYPE=INS;TANGLED=F;ES=>0>1;LV=0\tGT\t0\t1"
        ]
    );
    assert_eq!(
        native_record_lines("deletion_flubble.gfa", "HG1"),
        vec![
            "HG1#1#chr1\t1\t>0>1\tAG\tA\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>0>2,>0;VARTYPE=DEL;TANGLED=F;ES=>0>1;LV=0\tGT\t0\t1"
        ]
    );
}

#[test]
fn native_gfa_to_vcf_matches_lean_nested_leaf_fixtures() {
    assert_eq!(
        native_record_lines("nested_deletion.gfa", "HG1"),
        vec![
            "HG1#1#chr1\t2\t>1>4\tCT\tC\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1>3,>1;VARTYPE=DEL;TANGLED=F;ES=>1>4;LV=1\tGT\t0\t1\t."
        ]
    );
    assert_eq!(
        native_record_lines("nested_substitution_missing_outer.gfa", "HG1"),
        vec![
            "HG1#1#chr1\t3\t>1>4\tT\tG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>3,>6;VARTYPE=SUB;TANGLED=F;ES=>1>4;LV=1\tGT\t0\t1\t."
        ]
    );
}

#[test]
fn native_gfa_to_vcf_matches_lean_repeat_anchor_fixtures() {
    assert_eq!(
        native_record_lines("repeat_anchor_deletion.gfa", "HG1"),
        vec![
            "HG1#1#chr1\t1\t>0>2\tAA\tA\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>0>1,>0;VARTYPE=DEL;TANGLED=F;ES=>0>2;LV=0\tGT\t0\t1"
        ]
    );
    assert_eq!(
        native_record_lines("repeat_anchor_insertion.gfa", "HG1"),
        vec![
            "HG1#1#chr1\t1\t>0>2\tA\tAA\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>0,>0>1;VARTYPE=INS;TANGLED=F;ES=>0>2;LV=0\tGT\t0\t1"
        ]
    );
}

#[test]
fn native_gfa_to_vcf_matches_lean_complex_and_ordered_fixtures() {
    assert_eq!(
        native_record_lines("complex_substitution_span.gfa", "HG1"),
        vec![
            "HG1#1#chr1\t2\t>0>3\tCG\tTA\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1>2,>4>5;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0\tGT\t0\t1"
        ]
    );
    assert_eq!(
        native_record_lines("two_ordered_substitutions.gfa", "HG1"),
        vec![
            "HG1#1#chr1\t2\t>0>3\tC\tG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1,>2;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0\tGT\t0\t1",
            "HG1#1#chr1\t4\t>3>6\tA\tC\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>4,>5;VARTYPE=SUB;TANGLED=F;ES=>3>6;LV=0\tGT\t0\t1",
        ]
    );
}

#[test]
fn native_gfa_to_vcf_matches_lean_hairpin_and_linear_fixtures() {
    assert_eq!(
        native_record_lines("hairpin_inversion_subr.gfa", "ref"),
        vec![
            "ref\t2\t>1>5\tACGTA\tTACGT\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1>2>3>4>5,<5<4<3<2<1;VARTYPE=SUBR;TANGLED=F\tGT\t0\t1"
        ]
    );
    assert!(native_record_lines("linear_no_variant.gfa", "HG1").is_empty());
}

#[test]
fn native_gfa_to_vcf_rejects_unsupported_lean_failure_fixtures() {
    let temp = tempfile::tempdir().expect("tempdir");
    let ref_file = temp.path().join("refs.txt");
    fs::write(&ref_file, "HG1\n").expect("write reference file");

    for fixture in [
        "unsupported_overlap.gfa",
        "malformed_path_missing_overlaps.gfa",
    ] {
        let err = gfa_to_vcf_document(fixture_path(fixture), Some(&ref_file)).unwrap_err();
        assert!(
            matches!(err, Error::InvalidGfa { .. }),
            "expected InvalidGfa for {fixture}, got {err:?}"
        );
    }
}
