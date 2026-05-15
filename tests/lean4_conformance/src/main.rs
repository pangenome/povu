use std::collections::BTreeMap;
use std::env;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

#[derive(Clone, Debug)]
struct Fixture {
    id: &'static str,
    description: &'static str,
    gfa_path: PathBuf,
    reference_prefixes: Vec<&'static str>,
    expected: ExpectedOutcome,
}

#[derive(Clone, Debug)]
enum ExpectedOutcome {
    Vcf { check_record_order: bool },
    PovuFailure { reason: &'static str },
}

#[derive(Clone, Debug, Default)]
struct Options {
    repo_root: Option<PathBuf>,
    povu_bin: Option<PathBuf>,
    fixture_filter: Option<String>,
    skip_build: bool,
    skip_lean_build: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct NormalizedVcf {
    samples: Vec<String>,
    records: Vec<NormalizedRecord>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct NormalizedRecord {
    chrom: String,
    pos: u64,
    id: String,
    ref_allele: String,
    alt_alleles: Vec<String>,
    qual: String,
    filter: String,
    info: BTreeMap<String, String>,
    format: String,
    genotypes: BTreeMap<String, String>,
}

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let options = parse_args(env::args().skip(1))?;
    let repo_root = discover_repo_root(options.repo_root.as_deref())?;
    let fixtures = fixtures(&repo_root);
    let selected = select_fixtures(&fixtures, options.fixture_filter.as_deref())?;

    if !options.skip_lean_build {
        run_status(
            "lake build",
            Command::new("lake").arg("build").current_dir(&repo_root),
        )?;
    }

    let povu_bin = match options.povu_bin {
        Some(path) => absolutize(&repo_root, path),
        None => repo_root.join("bin").join("povu"),
    };

    if !options.skip_build {
        build_povu(&repo_root)?;
    } else if !povu_bin.exists() {
        return Err(format!(
            "--skip-build was passed, but povu binary does not exist: {}",
            povu_bin.display()
        ));
    }

    let mut failures = Vec::new();
    for fixture in selected {
        match run_fixture(&repo_root, &povu_bin, fixture) {
            Ok(()) => {
                println!("ok - {} ({})", fixture.id, fixture.description);
            }
            Err(err) => failures.push(err),
        }
    }

    if failures.is_empty() {
        println!("lean4 conformance: all fixtures passed");
        Ok(())
    } else {
        Err(failures.join("\n\n"))
    }
}

fn parse_args<I, S>(args: I) -> Result<Options, String>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let mut options = Options::default();
    let mut iter = args.into_iter().map(Into::into);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--help" | "-h" => {
                print_help();
                std::process::exit(0);
            }
            "--repo-root" => {
                options.repo_root = Some(PathBuf::from(next_value(&mut iter, "--repo-root")?));
            }
            "--povu-bin" => {
                options.povu_bin = Some(PathBuf::from(next_value(&mut iter, "--povu-bin")?));
            }
            "--fixture" => {
                options.fixture_filter = Some(next_value(&mut iter, "--fixture")?);
            }
            "--skip-build" => {
                options.skip_build = true;
            }
            "--skip-lean-build" => {
                options.skip_lean_build = true;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    Ok(options)
}

fn next_value<I>(iter: &mut I, flag: &str) -> Result<String, String>
where
    I: Iterator<Item = String>,
{
    iter.next()
        .ok_or_else(|| format!("missing value for {flag}"))
}

fn print_help() {
    println!(
        "povu Lean4 conformance harness\n\n\
         Usage:\n  \
         cargo run --manifest-path tests/lean4_conformance/Cargo.toml -- [options]\n\n\
         Options:\n  \
         --repo-root <path>       Repository root. Defaults to current directory or an ancestor.\n  \
         --povu-bin <path>        Existing povu binary to run. Defaults to <repo>/bin/povu.\n  \
         --fixture <id>           Run one fixture only.\n  \
         --skip-build             Do not configure/build the povu CLI before running fixtures.\n  \
         --skip-lean-build        Do not run lake build before running the Lean reference.\n"
    );
}

fn discover_repo_root(explicit: Option<&Path>) -> Result<PathBuf, String> {
    if let Some(root) = explicit {
        let root = env::current_dir()
            .map_err(|err| format!("failed to read current directory: {err}"))?
            .join(root);
        return validate_repo_root(root);
    }

    let mut dir =
        env::current_dir().map_err(|err| format!("failed to read current directory: {err}"))?;
    loop {
        if dir.join("lakefile.lean").exists() && dir.join("CMakeLists.txt").exists() {
            return validate_repo_root(dir);
        }
        if !dir.pop() {
            return Err(
                "could not find repository root containing lakefile.lean and CMakeLists.txt"
                    .to_string(),
            );
        }
    }
}

fn validate_repo_root(root: PathBuf) -> Result<PathBuf, String> {
    let root = root
        .canonicalize()
        .map_err(|err| format!("failed to canonicalize repo root {}: {err}", root.display()))?;
    if !root.join("lakefile.lean").exists() {
        return Err(format!("repo root lacks lakefile.lean: {}", root.display()));
    }
    if !root.join("CMakeLists.txt").exists() {
        return Err(format!(
            "repo root lacks CMakeLists.txt: {}",
            root.display()
        ));
    }
    Ok(root)
}

fn absolutize(repo_root: &Path, path: PathBuf) -> PathBuf {
    if path.is_absolute() {
        path
    } else {
        repo_root.join(path)
    }
}

fn fixtures(repo_root: &Path) -> Vec<Fixture> {
    let fixture_dir = repo_root.join("tests/lean4_conformance/fixtures");
    vec![
        Fixture {
            id: "minimal-substitution",
            description: "two haploid paths diverge at one base and rejoin",
            gfa_path: fixture_dir.join("minimal_substitution.gfa"),
            reference_prefixes: vec!["HG1"],
            expected: ExpectedOutcome::Vcf {
                check_record_order: false,
            },
        },
        Fixture {
            id: "insertion-flubble",
            description: "alternate path inserts one base between reference nodes",
            gfa_path: fixture_dir.join("insertion_flubble.gfa"),
            reference_prefixes: vec!["HG1"],
            expected: ExpectedOutcome::Vcf {
                check_record_order: false,
            },
        },
        Fixture {
            id: "deletion-flubble",
            description: "alternate path deletes one reference base with an anchor",
            gfa_path: fixture_dir.join("deletion_flubble.gfa"),
            reference_prefixes: vec!["HG1"],
            expected: ExpectedOutcome::Vcf {
                check_record_order: false,
            },
        },
        Fixture {
            id: "nested-deletion",
            description: "inner flubble emits a level-1 deletion and missing outer sample",
            gfa_path: fixture_dir.join("nested_deletion.gfa"),
            reference_prefixes: vec!["HG1"],
            expected: ExpectedOutcome::Vcf {
                check_record_order: false,
            },
        },
        Fixture {
            id: "hairpin-inversion-subr",
            description: "reverse traversal emits one SUBR hairpin inversion record",
            gfa_path: fixture_dir.join("hairpin_inversion_subr.gfa"),
            reference_prefixes: vec!["ref"],
            expected: ExpectedOutcome::Vcf {
                check_record_order: false,
            },
        },
        Fixture {
            id: "linear-no-variant",
            description: "identical haploid paths produce a header-only VCF",
            gfa_path: fixture_dir.join("linear_no_variant.gfa"),
            reference_prefixes: vec!["HG1"],
            expected: ExpectedOutcome::Vcf {
                check_record_order: false,
            },
        },
        Fixture {
            id: "two-ordered-substitutions",
            description: "two independent substitutions remain in deterministic VCF order",
            gfa_path: fixture_dir.join("two_ordered_substitutions.gfa"),
            reference_prefixes: vec!["HG1"],
            expected: ExpectedOutcome::Vcf {
                check_record_order: true,
            },
        },
        Fixture {
            id: "unsupported-overlap",
            description: "non-zero GFA link overlap is outside the accepted corpus subset",
            gfa_path: fixture_dir.join("unsupported_overlap.gfa"),
            reference_prefixes: vec!["HG1"],
            expected: ExpectedOutcome::PovuFailure {
                reason: "current gfa2vcf boundary rejects non-zero link overlaps",
            },
        },
        Fixture {
            id: "malformed-path-missing-overlaps",
            description: "malformed GFA path line missing its overlaps column is rejected",
            gfa_path: fixture_dir.join("malformed_path_missing_overlaps.gfa"),
            reference_prefixes: vec!["HG1"],
            expected: ExpectedOutcome::PovuFailure {
                reason: "current gfa2vcf boundary rejects malformed path records",
            },
        },
    ]
}

fn select_fixtures<'a>(
    fixtures: &'a [Fixture],
    fixture_filter: Option<&str>,
) -> Result<Vec<&'a Fixture>, String> {
    if let Some(id) = fixture_filter {
        let selected: Vec<_> = fixtures.iter().filter(|fixture| fixture.id == id).collect();
        if selected.is_empty() {
            return Err(format!("unknown fixture id: {id}"));
        }
        Ok(selected)
    } else {
        Ok(fixtures.iter().collect())
    }
}

fn build_povu(repo_root: &Path) -> Result<(), String> {
    let build_dir = repo_root.join("build").join("lean4-conformance");
    run_status(
        "cmake configure",
        Command::new("cmake")
            .arg("-S")
            .arg(repo_root)
            .arg("-B")
            .arg(&build_dir)
            .arg("-DPOVU_ENABLE_TESTING=ON")
            .current_dir(repo_root),
    )?;
    run_status(
        "cmake build povu",
        Command::new("cmake")
            .arg("--build")
            .arg(&build_dir)
            .arg("--target")
            .arg("povu")
            .current_dir(repo_root),
    )
}

fn run_status(label: &str, command: &mut Command) -> Result<(), String> {
    let status = command
        .status()
        .map_err(|err| format!("failed to run {label}: {err}"))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{label} failed with status {status}"))
    }
}

fn run_fixture(repo_root: &Path, povu_bin: &Path, fixture: &Fixture) -> Result<(), String> {
    let gfa = fs::read_to_string(&fixture.gfa_path).map_err(|err| {
        format!(
            "failed to read fixture {} ({}): {err}",
            fixture.id,
            fixture.gfa_path.display()
        )
    })?;

    let actual_output = run_povu(povu_bin, fixture)?;

    match &fixture.expected {
        ExpectedOutcome::Vcf { check_record_order } => run_vcf_fixture(
            repo_root,
            povu_bin,
            fixture,
            &gfa,
            actual_output,
            *check_record_order,
        ),
        ExpectedOutcome::PovuFailure { reason } => {
            run_expected_povu_failure(fixture, povu_bin, &gfa, actual_output, reason)
        }
    }
}

fn run_vcf_fixture(
    repo_root: &Path,
    povu_bin: &Path,
    fixture: &Fixture,
    gfa: &str,
    actual_output: Output,
    check_record_order: bool,
) -> Result<(), String> {
    let expected_output = run_lean_reference(repo_root, fixture)?;

    if !expected_output.status.success() {
        return Err(format_command_failure(
            fixture,
            "Lean reference",
            &lean_reference_command(repo_root, fixture),
            &expected_output,
        ));
    }
    if !actual_output.status.success() {
        return Err(format_command_failure(
            fixture,
            "povu gfa2vcf",
            &povu_command_display(povu_bin, fixture),
            &actual_output,
        ));
    }

    let expected_stdout = String::from_utf8_lossy(&expected_output.stdout);
    let actual_stdout = String::from_utf8_lossy(&actual_output.stdout);
    let expected = parse_vcf(&expected_stdout).map_err(|err| {
        format!(
            "Lean reference emitted invalid VCF for {}: {err}",
            fixture.id
        )
    })?;
    let actual = parse_vcf(&actual_stdout)
        .map_err(|err| format!("povu emitted invalid VCF for {}: {err}", fixture.id))?;

    let expected_for_compare = comparable_vcf(&expected, check_record_order);
    let actual_for_compare = comparable_vcf(&actual, check_record_order);

    if expected_for_compare == actual_for_compare {
        Ok(())
    } else {
        Err(render_mismatch(
            fixture,
            &povu_command_display(povu_bin, fixture),
            &gfa,
            &expected,
            &actual,
            &actual_stdout,
            &String::from_utf8_lossy(&actual_output.stderr),
        ))
    }
}

fn run_expected_povu_failure(
    fixture: &Fixture,
    povu_bin: &Path,
    gfa: &str,
    actual_output: Output,
    reason: &str,
) -> Result<(), String> {
    if !actual_output.status.success() {
        return actual_output.status.code().map(|_| ()).ok_or_else(|| {
            render_uncontrolled_failure(
                fixture,
                &povu_command_display(povu_bin, fixture),
                gfa,
                &String::from_utf8_lossy(&actual_output.stdout),
                &String::from_utf8_lossy(&actual_output.stderr),
                reason,
                &actual_output,
            )
        });
    }

    Err(render_unexpected_success(
        fixture,
        &povu_command_display(povu_bin, fixture),
        gfa,
        &String::from_utf8_lossy(&actual_output.stdout),
        &String::from_utf8_lossy(&actual_output.stderr),
        reason,
    ))
}

fn run_lean_reference(repo_root: &Path, fixture: &Fixture) -> Result<Output, String> {
    Command::new("lake")
        .arg("env")
        .arg("lean")
        .arg("--run")
        .arg(repo_root.join("tests/lean4_conformance/lean_reference.lean"))
        .arg(fixture.id)
        .current_dir(repo_root)
        .output()
        .map_err(|err| format!("failed to run Lean reference for {}: {err}", fixture.id))
}

fn run_povu(povu_bin: &Path, fixture: &Fixture) -> Result<Output, String> {
    let mut command = Command::new(povu_bin);
    command
        .arg("gfa2vcf")
        .arg("-i")
        .arg(&fixture.gfa_path)
        .arg("-t")
        .arg("1");
    for prefix in &fixture.reference_prefixes {
        command.arg("-P").arg(prefix);
    }
    command
        .output()
        .map_err(|err| format!("failed to run povu for {}: {err}", fixture.id))
}

fn lean_reference_command(repo_root: &Path, fixture: &Fixture) -> Vec<String> {
    vec![
        "lake".to_string(),
        "env".to_string(),
        "lean".to_string(),
        "--run".to_string(),
        repo_root
            .join("tests/lean4_conformance/lean_reference.lean")
            .display()
            .to_string(),
        fixture.id.to_string(),
    ]
}

fn povu_command_display(povu_bin: &Path, fixture: &Fixture) -> Vec<String> {
    let mut command = vec![
        povu_bin.display().to_string(),
        "gfa2vcf".to_string(),
        "-i".to_string(),
        fixture.gfa_path.display().to_string(),
        "-t".to_string(),
        "1".to_string(),
    ];
    for prefix in &fixture.reference_prefixes {
        command.push("-P".to_string());
        command.push((*prefix).to_string());
    }
    command
}

fn format_command_failure(
    fixture: &Fixture,
    label: &str,
    command: &[String],
    output: &Output,
) -> String {
    format!(
        "{label} failed for fixture {id}\ncommand: {command}\nstatus: {status}\nstdout:\n{stdout}\nstderr:\n{stderr}",
        id = fixture.id,
        command = shell_join(command),
        status = output.status,
        stdout = String::from_utf8_lossy(&output.stdout),
        stderr = String::from_utf8_lossy(&output.stderr)
    )
}

fn parse_vcf(text: &str) -> Result<NormalizedVcf, String> {
    let mut samples = Vec::new();
    let mut records = Vec::new();
    let mut saw_header = false;

    for (line_idx, line) in text.lines().enumerate() {
        let line_no = line_idx + 1;
        if line.trim().is_empty() || line.starts_with("##") {
            continue;
        }

        if line.starts_with("#CHROM") {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() < 9 {
                return Err(format!(
                    "line {line_no}: VCF column header has fewer than 9 fields"
                ));
            }
            samples = fields[9..]
                .iter()
                .map(|field| (*field).to_string())
                .collect();
            saw_header = true;
            continue;
        }

        if !saw_header {
            return Err(format!(
                "line {line_no}: record appeared before #CHROM header"
            ));
        }

        let fields: Vec<_> = line.split('\t').collect();
        let expected_fields = 9 + samples.len();
        if fields.len() != expected_fields {
            return Err(format!(
                "line {line_no}: expected {expected_fields} tab-separated fields, found {}",
                fields.len()
            ));
        }

        let pos = fields[1]
            .parse::<u64>()
            .map_err(|err| format!("line {line_no}: invalid POS '{}': {err}", fields[1]))?;
        let genotypes = samples
            .iter()
            .cloned()
            .zip(fields[9..].iter().map(|field| (*field).to_string()))
            .collect();

        records.push(NormalizedRecord {
            chrom: fields[0].to_string(),
            pos,
            id: fields[2].to_string(),
            ref_allele: fields[3].to_string(),
            alt_alleles: split_list(fields[4]),
            qual: fields[5].to_string(),
            filter: fields[6].to_string(),
            info: parse_info(fields[7], line_no)?,
            format: fields[8].to_string(),
            genotypes,
        });
    }

    if !saw_header {
        return Err("missing #CHROM header".to_string());
    }

    Ok(NormalizedVcf { samples, records })
}

fn comparable_vcf(vcf: &NormalizedVcf, check_record_order: bool) -> NormalizedVcf {
    let mut comparable = vcf.clone();
    if !check_record_order {
        sort_records(&mut comparable.records);
    }
    comparable
}

fn sort_records(records: &mut [NormalizedRecord]) {
    records.sort_by(|left, right| {
        left.chrom
            .cmp(&right.chrom)
            .then(left.pos.cmp(&right.pos))
            .then(left.id.cmp(&right.id))
            .then(left.ref_allele.cmp(&right.ref_allele))
            .then(left.alt_alleles.cmp(&right.alt_alleles))
    });
}

fn split_list(value: &str) -> Vec<String> {
    if value == "." || value.is_empty() {
        Vec::new()
    } else {
        value.split(',').map(str::to_string).collect()
    }
}

fn parse_info(value: &str, line_no: usize) -> Result<BTreeMap<String, String>, String> {
    let mut info = BTreeMap::new();
    if value == "." || value.is_empty() {
        return Ok(info);
    }

    for entry in value.split(';') {
        if entry.is_empty() {
            return Err(format!("line {line_no}: empty INFO entry"));
        }
        let (key, field_value) = entry
            .split_once('=')
            .map_or((entry, ""), |(key, field_value)| (key, field_value));
        if key.is_empty() {
            return Err(format!("line {line_no}: INFO entry has empty key"));
        }
        if info
            .insert(key.to_string(), field_value.to_string())
            .is_some()
        {
            return Err(format!("line {line_no}: duplicate INFO key '{key}'"));
        }
    }
    Ok(info)
}

fn render_mismatch(
    fixture: &Fixture,
    command: &[String],
    gfa: &str,
    expected: &NormalizedVcf,
    actual: &NormalizedVcf,
    raw_stdout: &str,
    raw_stderr: &str,
) -> String {
    let mut message = String::new();
    writeln!(
        &mut message,
        "VCF semantic mismatch for fixture '{}' ({})",
        fixture.id, fixture.description
    )
    .unwrap();
    writeln!(&mut message, "command: {}", shell_join(command)).unwrap();
    writeln!(
        &mut message,
        "input fixture: {}",
        fixture.gfa_path.display()
    )
    .unwrap();
    writeln!(&mut message, "\n--- input GFA ---\n{gfa}").unwrap();
    writeln!(
        &mut message,
        "--- expected semantic VCF (Lean) ---\n{}",
        format_vcf(expected)
    )
    .unwrap();
    writeln!(
        &mut message,
        "--- actual semantic VCF (povu) ---\n{}",
        format_vcf(actual)
    )
    .unwrap();
    writeln!(&mut message, "--- raw povu stdout ---\n{raw_stdout}").unwrap();
    if !raw_stderr.trim().is_empty() {
        writeln!(&mut message, "--- raw povu stderr ---\n{raw_stderr}").unwrap();
    }
    message
}

fn render_unexpected_success(
    fixture: &Fixture,
    command: &[String],
    gfa: &str,
    raw_stdout: &str,
    raw_stderr: &str,
    reason: &str,
) -> String {
    let mut message = String::new();
    writeln!(
        &mut message,
        "povu succeeded for fixture '{}' but a failure was expected",
        fixture.id
    )
    .unwrap();
    writeln!(&mut message, "reason: {reason}").unwrap();
    writeln!(&mut message, "command: {}", shell_join(command)).unwrap();
    writeln!(
        &mut message,
        "input fixture: {}",
        fixture.gfa_path.display()
    )
    .unwrap();
    writeln!(&mut message, "\n--- input GFA ---\n{gfa}").unwrap();
    writeln!(&mut message, "--- raw povu stdout ---\n{raw_stdout}").unwrap();
    if !raw_stderr.trim().is_empty() {
        writeln!(&mut message, "--- raw povu stderr ---\n{raw_stderr}").unwrap();
    }
    message
}

fn render_uncontrolled_failure(
    fixture: &Fixture,
    command: &[String],
    gfa: &str,
    raw_stdout: &str,
    raw_stderr: &str,
    reason: &str,
    output: &Output,
) -> String {
    let mut message = String::new();
    writeln!(
        &mut message,
        "povu failed without a controlled exit code for expected-failure fixture '{}'",
        fixture.id
    )
    .unwrap();
    writeln!(&mut message, "reason: {reason}").unwrap();
    writeln!(&mut message, "command: {}", shell_join(command)).unwrap();
    writeln!(&mut message, "status: {}", output.status).unwrap();
    writeln!(
        &mut message,
        "input fixture: {}",
        fixture.gfa_path.display()
    )
    .unwrap();
    writeln!(&mut message, "\n--- input GFA ---\n{gfa}").unwrap();
    writeln!(&mut message, "--- raw povu stdout ---\n{raw_stdout}").unwrap();
    if !raw_stderr.trim().is_empty() {
        writeln!(&mut message, "--- raw povu stderr ---\n{raw_stderr}").unwrap();
    }
    message
}

fn format_vcf(vcf: &NormalizedVcf) -> String {
    let mut out = String::new();
    writeln!(&mut out, "samples: {}", vcf.samples.join(", ")).unwrap();
    if vcf.records.is_empty() {
        out.push_str("records: <none>\n");
        return out;
    }

    out.push_str("records:\n");
    for record in &vcf.records {
        writeln!(
            &mut out,
            "  {}:{} {} {} -> {}",
            record.chrom,
            record.pos,
            record.id,
            record.ref_allele,
            if record.alt_alleles.is_empty() {
                ".".to_string()
            } else {
                record.alt_alleles.join(",")
            }
        )
        .unwrap();
        writeln!(
            &mut out,
            "    qual/filter/format: {}/{}/{}",
            record.qual, record.filter, record.format
        )
        .unwrap();
        out.push_str("    info:\n");
        for (key, value) in &record.info {
            writeln!(&mut out, "      {key}={value}").unwrap();
        }
        out.push_str("    genotypes:\n");
        for (sample, genotype) in &record.genotypes {
            writeln!(&mut out, "      {sample}={genotype}").unwrap();
        }
    }
    out
}

fn shell_join(parts: &[String]) -> String {
    parts
        .iter()
        .map(|part| shell_quote(part))
        .collect::<Vec<_>>()
        .join(" ")
}

fn shell_quote(part: &str) -> String {
    if part.is_empty()
        || part
            .chars()
            .any(|ch| ch.is_whitespace() || matches!(ch, '\'' | '"' | '$' | '\\' | ';' | '&' | '|'))
    {
        format!("'{}'", part.replace('\'', "'\\''"))
    } else {
        part.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\n";

    #[test]
    fn normalizes_info_field_order() {
        let left = format!(
            "{}HG1#1#chr1\t2\t>0>3\tC\tG\t60\tPASS\tAC=1;AF=0.5;AN=2;NS=2;AT=>1,>2;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0\tGT\t0\t1\n",
            HEADER
        );
        let right = format!(
            "{}HG1#1#chr1\t2\t>0>3\tC\tG\t60\tPASS\tLV=0;ES=>0>3;TANGLED=F;VARTYPE=SUB;AT=>1,>2;NS=2;AN=2;AF=0.5;AC=1\tGT\t0\t1\n",
            HEADER
        );
        assert_eq!(parse_vcf(&left).unwrap(), parse_vcf(&right).unwrap());
    }

    #[test]
    fn mismatch_diagnostic_includes_fixture_and_semantic_outputs() {
        let fixture = Fixture {
            id: "unit",
            description: "diagnostic coverage",
            gfa_path: PathBuf::from("tests/lean4_conformance/fixtures/unit.gfa"),
            reference_prefixes: vec!["HG1"],
            expected: ExpectedOutcome::Vcf {
                check_record_order: false,
            },
        };
        let expected = parse_vcf(&format!(
            "{}HG1#1#chr1\t2\t>0>3\tC\tG\t60\tPASS\tAC=1\tGT\t0\t1\n",
            HEADER
        ))
        .unwrap();
        let actual = parse_vcf(&format!(
            "{}HG1#1#chr1\t2\t>0>3\tC\tT\t60\tPASS\tAC=1\tGT\t0\t1\n",
            HEADER
        ))
        .unwrap();
        let diagnostic = render_mismatch(
            &fixture,
            &["povu".to_string(), "gfa2vcf".to_string()],
            "S\t0\tA\n",
            &expected,
            &actual,
            "raw stdout",
            "",
        );

        assert!(diagnostic.contains("input fixture: tests/lean4_conformance/fixtures/unit.gfa"));
        assert!(diagnostic.contains("--- input GFA ---"));
        assert!(diagnostic.contains("--- expected semantic VCF (Lean) ---"));
        assert!(diagnostic.contains("--- actual semantic VCF (povu) ---"));
        assert!(diagnostic.contains("C -> G"));
        assert!(diagnostic.contains("C -> T"));
    }

    #[test]
    fn comparable_vcf_can_preserve_or_normalize_record_order() {
        let first = "HG1#1#chr1\t4\t>3>6\tA\tC\t60\tPASS\tAC=1\tGT\t0\t1\n";
        let second = "HG1#1#chr1\t2\t>0>3\tC\tG\t60\tPASS\tAC=1\tGT\t0\t1\n";
        let parsed = parse_vcf(&format!("{HEADER}{first}{second}")).unwrap();

        assert_eq!(parsed.records[0].pos, 4);
        assert_eq!(comparable_vcf(&parsed, true).records[0].pos, 4);
        assert_eq!(comparable_vcf(&parsed, false).records[0].pos, 2);
    }
}
