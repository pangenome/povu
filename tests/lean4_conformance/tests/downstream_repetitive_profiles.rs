use std::path::Path;
use std::process::Command;

use povu_lean4_conformance::downstream_profiles::{
    run_downstream_repetitive_corpus, run_fixture_profile, ProfileOutput,
};

fn repo_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("lean4 conformance crate lives under tests/lean4_conformance")
}

fn vcf_text(fixture_id: &str, profile_name: &str) -> String {
    match run_fixture_profile(repo_root(), fixture_id, profile_name).unwrap() {
        ProfileOutput::Vcf { vcf, .. } => vcf,
        ProfileOutput::Reject { reason } => {
            panic!("{fixture_id}/{profile_name} rejected unexpectedly: {reason}")
        }
    }
}

fn provenance_text(fixture_id: &str, profile_name: &str) -> String {
    match run_fixture_profile(repo_root(), fixture_id, profile_name).unwrap() {
        ProfileOutput::Vcf {
            provenance: Some(provenance),
            ..
        } => provenance,
        ProfileOutput::Vcf {
            provenance: None, ..
        } => panic!("{fixture_id}/{profile_name} did not emit provenance"),
        ProfileOutput::Reject { reason } => {
            panic!("{fixture_id}/{profile_name} rejected unexpectedly: {reason}")
        }
    }
}

#[test]
fn tandem_repeat_left_normalization_shifts_indels_and_emits_provenance() {
    let raw = vcf_text("tandem-repeat-left-normalization", "raw-graph");
    let normalized = vcf_text("tandem-repeat-left-normalization", "left-normalized");
    let provenance = provenance_text("tandem-repeat-left-normalization", "left-normalized");

    assert!(raw.contains("\t4\t>3>4\tA\tAA\t"));
    assert!(raw.contains("\t4\t>3>5\tAA\tA\t"));
    assert!(normalized.contains("\t1\t>3>4:norm\tA\tAA\t"));
    assert!(normalized.contains("\t1\t>3>5:norm\tAA\tA\t"));
    assert!(normalized.contains("LEFT_NORMALIZED=T"));
    assert!(provenance.contains("\"operation\":\"coordinate-shift\""));
    assert!(provenance.contains("\"raw_pos\":4"));
    assert!(provenance.contains("\"output_pos\":1"));
}

#[test]
fn popped_profile_rescues_child_of_popped_parent() {
    let top_level = vcf_text("popped-parent-child-rescue", "top-level-only");
    let popped = vcf_text("popped-parent-child-rescue", "popped");
    let provenance = provenance_text("popped-parent-child-rescue", "popped");

    assert!(top_level.contains("\t1\t>0>5:top\tAAAAA\tA\t"));
    assert!(!top_level.contains(">2>4:rescued"));
    assert!(!popped.contains("\t1\t>0>5\tAAAAA\tA\t"));
    assert!(popped.contains("\t3\t>2>4:rescued\tA\tG\t"));
    assert!(popped.contains("RESCUED_CHILD=T"));
    assert!(provenance.contains("\"operation\":\"pop-parent\""));
    assert!(provenance.contains("\"operation\":\"rescue-child\""));
}

#[test]
fn decomposed_profile_splits_short_alt_and_passes_long_alt() {
    let decomposed = vcf_text("vcfwave-complex-decomposition", "decomposed");
    let provenance = provenance_text("vcfwave-complex-decomposition", "decomposed");

    assert!(decomposed.contains("\t11\t>9>14:1:snp1\tC\tT\t"));
    assert!(decomposed.contains("\t13\t>9>14:1:snp2\tT\tA\t"));
    assert!(decomposed.contains("\t10\t>9>14:2:passthrough\tACGT\tACGTACTACGTACGTA\t"));
    assert!(decomposed.contains("DECOMPOSED=T"));
    assert!(decomposed.contains("PASSTHROUGH=T"));
    assert!(decomposed.contains("\tGT\t0\t1\t.\n"));
    assert!(decomposed.contains("\tGT\t0\t.\t1\n"));
    assert!(provenance.contains("\"operation\":\"one-to-many-decomposition\""));
    assert!(provenance.contains("\">9>14:1:snp1\""));
    assert!(provenance.contains("\">9>14:1:snp2\""));
}

#[test]
fn nested_child_inside_insertion_is_expected_reject() {
    let raw = run_fixture_profile(repo_root(), "nested-child-inside-insertion", "raw-graph")
        .expect("raw profile should produce a controlled outcome");
    let inserted = run_fixture_profile(
        repo_root(),
        "nested-child-inside-insertion",
        "inserted-path-raw",
    )
    .expect("inserted path profile should produce a controlled outcome");

    match raw {
        ProfileOutput::Reject { reason } => {
            assert!(reason.contains("no selected-reference coordinate"))
        }
        ProfileOutput::Vcf { .. } => panic!("raw profile must reject this nested insertion child"),
    }
    match inserted {
        ProfileOutput::Reject { reason } => {
            assert!(reason.contains("inserted-path coordinate profile"))
        }
        ProfileOutput::Vcf { .. } => {
            panic!("inserted-path profile must reject until coordinates are specified")
        }
    }
}

#[test]
fn subr_inversion_decomposition_preserves_graph_native_record() {
    let raw = vcf_text("subr-inversion-preservation", "raw-graph");
    let decomposed = vcf_text("subr-inversion-preservation", "decomposed");
    let provenance = provenance_text("subr-inversion-preservation", "decomposed");

    assert!(raw.contains("VARTYPE=SUBR"));
    assert!(decomposed.contains("VARTYPE=SUBR"));
    assert!(decomposed.contains("SUBR_ORIGIN=T"));
    assert!(!decomposed.contains("<INV>"));
    assert!(!decomposed.contains(";ES="));
    assert!(!decomposed.contains(";LV="));
    assert!(provenance.contains("\"SUBR_PRESERVED\""));
}

#[test]
fn downstream_repetitive_manifest_profiles_pass() {
    run_downstream_repetitive_corpus(repo_root(), None).unwrap();
}

#[test]
fn cli_runs_downstream_repetitive_manifest_filter() {
    let output = Command::new(env!("CARGO_BIN_EXE_povu-lean4-conformance"))
        .arg("--repo-root")
        .arg(repo_root())
        .arg("--downstream-repetitive")
        .arg("--fixture")
        .arg("tandem-repeat-left-normalization")
        .output()
        .expect("failed to run conformance CLI");

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}
