use povu::{detect_flubble_stack, gfa_to_vcf, gfa_to_vcf_document, Error, FlubbleCandidate};
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
