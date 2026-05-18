use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use serde_json::Value;

#[test]
fn downstream_repetitive_fixture_assets_are_complete() {
    validate_downstream_repetitive_assets().unwrap();
}

fn validate_downstream_repetitive_assets() -> Result<(), String> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("fixtures/downstream_repetitive");
    let manifest_path = root.join("manifest.json");
    let manifest = read_json(&manifest_path)?;
    expect_string(&manifest, "schema", &manifest_path)?
        .eq("povu.downstream_repetitive.fixtures.v1")
        .then_some(())
        .ok_or_else(|| format!("{} has unexpected schema", manifest_path.display()))?;

    let fixtures = expect_array(&manifest, "fixtures", &manifest_path)?;
    if fixtures.is_empty() {
        return Err(format!(
            "{} does not list any fixtures",
            manifest_path.display()
        ));
    }

    let mut coverage = BTreeSet::new();
    let mut saw_coordinate_shift_provenance = false;
    let mut saw_one_to_many_provenance = false;
    let mut saw_expected_reject = false;

    for fixture in fixtures {
        let fixture_id = expect_string(fixture, "id", &manifest_path)?;
        for tag in expect_array(fixture, "class_tags", &manifest_path)? {
            coverage.insert(expect_value_string(tag, "class tag", &manifest_path)?.to_string());
        }

        let input_gfa = root.join(expect_string(fixture, "input_gfa", &manifest_path)?);
        require_file(&input_gfa)?;
        if let Some(reference_fasta) = optional_string(fixture, "reference_fasta") {
            require_file(&root.join(reference_fasta))?;
        }

        let profiles = expect_array(fixture, "profiles", &manifest_path)?;
        let mut raw_artifacts = BTreeSet::new();
        let mut downstream_artifacts = BTreeSet::new();
        let mut has_raw = false;
        let mut has_downstream = false;

        for profile in profiles {
            let profile_name = expect_string(profile, "name", &manifest_path)?;
            let profile_kind = expect_string(profile, "kind", &manifest_path)?;
            let expected = expect_string(profile, "expected", &manifest_path)?;

            let artifact = match expected {
                "vcf" => {
                    let rel = expect_string(profile, "path", &manifest_path)?;
                    let path = root.join(rel);
                    validate_vcf(&path, fixture_id, profile_name)?;

                    if let Some(provenance_rel) = optional_string(profile, "provenance") {
                        let stats = validate_provenance(&root.join(provenance_rel), fixture_id)?;
                        match optional_string(profile, "rewrite") {
                            Some("coordinate-shift") => {
                                if !stats.saw_coordinate_shift {
                                    return Err(format!(
                                        "{fixture_id}/{profile_name} declares coordinate-shift but provenance has no raw/output POS shift"
                                    ));
                                }
                                saw_coordinate_shift_provenance = true;
                            }
                            Some("one-to-many") => {
                                if !stats.saw_one_to_many {
                                    return Err(format!(
                                        "{fixture_id}/{profile_name} declares one-to-many but provenance has no multi-output mapping"
                                    ));
                                }
                                saw_one_to_many_provenance = true;
                            }
                            Some("pass-through") | None => {}
                            Some(other) => {
                                return Err(format!(
                                    "{fixture_id}/{profile_name} uses unknown rewrite marker {other}"
                                ));
                            }
                        }
                    }

                    path
                }
                "reject" => {
                    let rel = expect_string(profile, "reject_path", &manifest_path)?;
                    let path = root.join(rel);
                    validate_expected_reject(&path, fixture_id, profile_name)?;
                    saw_expected_reject = true;
                    path
                }
                other => {
                    return Err(format!(
                        "{fixture_id}/{profile_name} has unsupported expected kind {other}"
                    ));
                }
            };

            if profile_kind == "raw" {
                has_raw = true;
                raw_artifacts.insert(artifact);
            } else {
                has_downstream = true;
                downstream_artifacts.insert(artifact);
            }
        }

        if !has_raw {
            return Err(format!("{fixture_id} has no raw profile expectation"));
        }
        if !has_downstream {
            return Err(format!(
                "{fixture_id} has no downstream profile expectation distinct from raw output"
            ));
        }
        for raw_artifact in raw_artifacts {
            if downstream_artifacts.contains(&raw_artifact) {
                return Err(format!(
                    "{fixture_id} reuses raw artifact {} as downstream output",
                    raw_artifact.display()
                ));
            }
        }
    }

    for required in [
        "tandem-repeat-normalization",
        "parent-popping-child-rescue",
        "complex-decomposition-pass-through",
        "nested-child-inside-insertion",
        "subr-inversion-preservation",
    ] {
        if !coverage.contains(required) {
            return Err(format!(
                "downstream repetitive corpus lacks coverage tag {required}"
            ));
        }
    }

    if !saw_coordinate_shift_provenance {
        return Err("coordinate-shifting rewrite provenance was not validated".to_string());
    }
    if !saw_one_to_many_provenance {
        return Err("one-to-many rewrite provenance was not validated".to_string());
    }
    if !saw_expected_reject {
        return Err("expected-reject policy assets were not validated".to_string());
    }

    Ok(())
}

#[derive(Default)]
struct ProvenanceStats {
    saw_coordinate_shift: bool,
    saw_one_to_many: bool,
}

fn validate_vcf(path: &Path, fixture_id: &str, profile_name: &str) -> Result<(), String> {
    require_file(path)?;
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    let mut saw_column_header = false;
    let mut header_width = 0usize;
    let mut record_count = 0usize;

    for (idx, line) in text.lines().enumerate() {
        let line_no = idx + 1;
        if line.trim().is_empty() || line.starts_with("##") {
            continue;
        }
        if line.starts_with("#CHROM") {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() < 10 {
                return Err(format!(
                    "{}:{line_no}: VCF header must include FORMAT and at least one sample",
                    path.display()
                ));
            }
            saw_column_header = true;
            header_width = fields.len();
            continue;
        }

        if !saw_column_header {
            return Err(format!(
                "{}:{line_no}: data record appears before #CHROM header",
                path.display()
            ));
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != header_width {
            return Err(format!(
                "{}:{line_no}: expected {header_width} VCF columns, found {}",
                path.display(),
                fields.len()
            ));
        }
        if fields[0].is_empty()
            || fields[2].is_empty()
            || fields[3].is_empty()
            || fields[4].is_empty()
        {
            return Err(format!(
                "{}:{line_no}: CHROM, ID, REF, and ALT must be nonempty",
                path.display()
            ));
        }
        let pos: u64 = fields[1].parse().map_err(|err| {
            format!(
                "{}:{line_no}: POS is not a positive integer: {err}",
                path.display()
            )
        })?;
        if pos == 0 {
            return Err(format!(
                "{}:{line_no}: POS must be positive",
                path.display()
            ));
        }
        if fields[5] != "60" || fields[6] != "PASS" || fields[8] != "GT" {
            return Err(format!(
                "{}:{line_no}: expected QUAL=60, FILTER=PASS, FORMAT=GT",
                path.display()
            ));
        }

        let info = parse_info(fields[7]);
        for required in ["AC", "AF", "AN", "NS", "AT", "VARTYPE", "TANGLED"] {
            if !info.contains_key(required) {
                return Err(format!(
                    "{}:{line_no}: INFO lacks required field {required}",
                    path.display()
                ));
            }
        }
        if fields[4].contains("<INV>") {
            return Err(format!(
                "{}:{line_no}: SUBR preservation fixtures must not use symbolic <INV>",
                path.display()
            ));
        }
        if info.get("VARTYPE").is_some_and(|value| *value == "SUBR")
            && (info.contains_key("ES") || info.contains_key("LV"))
        {
            return Err(format!(
                "{}:{line_no}: SUBR records must not carry flubble ES/LV fields",
                path.display()
            ));
        }
        record_count += 1;
    }

    if !saw_column_header {
        return Err(format!("{} lacks #CHROM header", path.display()));
    }
    if record_count == 0 {
        return Err(format!(
            "{fixture_id}/{profile_name} VCF {} contains no records",
            path.display()
        ));
    }

    Ok(())
}

fn parse_info(text: &str) -> BTreeMap<&str, &str> {
    text.split(';')
        .filter(|field| !field.is_empty())
        .map(|field| field.split_once('=').unwrap_or((field, "T")))
        .collect()
}

fn validate_expected_reject(
    path: &Path,
    fixture_id: &str,
    profile_name: &str,
) -> Result<(), String> {
    require_file(path)?;
    let value = read_json(path)?;
    if expect_string(&value, "fixture_id", path)? != fixture_id {
        return Err(format!("{} has wrong fixture_id", path.display()));
    }
    if expect_string(&value, "profile", path)? != profile_name {
        return Err(format!("{} has wrong profile", path.display()));
    }
    if expect_string(&value, "status", path)? != "expected-reject" {
        return Err(format!("{} is not marked expected-reject", path.display()));
    }
    if expect_string(&value, "reason", path)?.is_empty() {
        return Err(format!("{} has an empty reject reason", path.display()));
    }
    if expect_string(&value, "required_behavior", path)?.is_empty() {
        return Err(format!("{} has empty required_behavior", path.display()));
    }
    Ok(())
}

fn validate_provenance(path: &Path, fixture_id: &str) -> Result<ProvenanceStats, String> {
    require_file(path)?;
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    let mut stats = ProvenanceStats::default();
    let mut lines = 0usize;

    for (idx, line) in text.lines().enumerate() {
        let line_no = idx + 1;
        if line.trim().is_empty() {
            continue;
        }
        lines += 1;
        let value: Value = serde_json::from_str(line).map_err(|err| {
            format!(
                "{}:{line_no}: invalid provenance JSONL record: {err}",
                path.display()
            )
        })?;
        if expect_string(&value, "fixture_id", path)? != fixture_id {
            return Err(format!("{}:{line_no}: wrong fixture_id", path.display()));
        }
        for key in [
            "profile",
            "operation",
            "raw_record_id",
            "raw_chrom",
            "raw_ref",
            "raw_alt",
        ] {
            if expect_string(&value, key, path)?.is_empty() {
                return Err(format!(
                    "{}:{line_no}: provenance field {key} must be nonempty",
                    path.display()
                ));
            }
        }
        expect_u64(&value, "raw_alt_index", path)?;
        let raw_pos = expect_u64(&value, "raw_pos", path)?;
        expect_array(&value, "raw_at", path)?;
        let outputs = expect_array(&value, "output_record_ids", path)?;
        let flags = expect_array(&value, "flags", path)?;
        expect_object(&value, "profile_config", path)?;

        if outputs.len() > 1 {
            stats.saw_one_to_many = true;
        }
        if let Some(output_pos) = optional_u64(&value, "output_pos") {
            if output_pos != raw_pos {
                stats.saw_coordinate_shift = true;
            }
        }
        for flag in flags {
            expect_value_string(flag, "flag", path)?;
        }
    }

    if lines == 0 {
        return Err(format!("{} has no provenance records", path.display()));
    }

    Ok(stats)
}

fn require_file(path: &Path) -> Result<(), String> {
    if path.is_file() {
        Ok(())
    } else {
        Err(format!(
            "required fixture asset is missing: {}",
            path.display()
        ))
    }
}

fn read_json(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    serde_json::from_str(&text).map_err(|err| format!("invalid JSON in {}: {err}", path.display()))
}

fn expect_string<'a>(value: &'a Value, key: &str, path: &Path) -> Result<&'a str, String> {
    value
        .get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("{} lacks string field {key}", path.display()))
}

fn optional_string<'a>(value: &'a Value, key: &str) -> Option<&'a str> {
    value.get(key).and_then(Value::as_str)
}

fn expect_value_string<'a>(value: &'a Value, label: &str, path: &Path) -> Result<&'a str, String> {
    value
        .as_str()
        .ok_or_else(|| format!("{} has non-string {label}", path.display()))
}

fn expect_array<'a>(value: &'a Value, key: &str, path: &Path) -> Result<&'a [Value], String> {
    value
        .get(key)
        .and_then(Value::as_array)
        .map(Vec::as_slice)
        .ok_or_else(|| format!("{} lacks array field {key}", path.display()))
}

fn expect_object<'a>(
    value: &'a Value,
    key: &str,
    path: &Path,
) -> Result<&'a serde_json::Map<String, Value>, String> {
    value
        .get(key)
        .and_then(Value::as_object)
        .ok_or_else(|| format!("{} lacks object field {key}", path.display()))
}

fn expect_u64(value: &Value, key: &str, path: &Path) -> Result<u64, String> {
    value
        .get(key)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("{} lacks unsigned integer field {key}", path.display()))
}

fn optional_u64(value: &Value, key: &str) -> Option<u64> {
    value.get(key).and_then(Value::as_u64)
}
