import csv
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent
PHASES = [
    ("original_default", pathlib.Path("/Users/vinky/research-comps-clean/production_checks/numerical_stability"), 0.005, 0.001, True),
    ("tight_tolerance", pathlib.Path("/tmp/psmelding-tight-tolerance-20261022"), 0.0001, 0.001, True),
    ("fixed_backend", pathlib.Path("/tmp/psmelding-fixed-backend-20261022"), 0.0001, 0.001, False),
    ("tight_inner", pathlib.Path("/tmp/psmelding-tight-inner-20261022"), 0.0001, 0.00001, False),
]

def read_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))

def write_rows(name, rows):
    if rows:
        with (ROOT / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

stability, attempts, phases = [], [], []
for phase, path, outer, inner, internal in PHASES:
    report = path / "reports" / "stability_by_prior.csv"
    if not report.exists():
        continue
    metadata = dict(phase=phase, outer_tolerance=outer, inner_tolerance=inner,
                    internal_opt=internal)
    for row in read_rows(report):
        stability.append(dict(metadata, **row))
    for row in read_rows(path / "reports" / "all_attempts.csv"):
        attempts.append(dict(metadata, **row))
    phases.append(dict(metadata, execution_or_reference_root=str(path)))

write_rows("phase_stability.csv", stability)
write_rows("all_phase_attempts.csv", attempts)
write_rows("phases.csv", phases)
score_ranges = []
for phase in phases:
    for dataset in ("baseline", "clustering_delta2"):
        selected = [row for row in attempts
                    if row['phase'] == phase['phase'] and row['dataset'] == dataset]
        result = dict(phase=phase['phase'], dataset=dataset)
        for metric in ("MSE", "CRPS_Z_joint"):
            values = [float(row[metric]) for row in selected]
            result[metric + '_min'] = min(values)
            result[metric + '_max'] = max(values)
        score_ranges.append(result)
write_rows("score_ranges.csv", score_ranges)
print(f"Summarized {len(phases)} phases, {len(attempts)} retained attempts")

if len(phases) == 4:
    lines = [
        "NUMERICAL TOLERANCE FOLLOW-UP -- 2026-10-22",
        "",
        "Study launch remains pending; the main seven R scripts were not changed",
        "18 new controlled fits are retained, alongside the six original reference fits",
        "No original fits were replaced and no full simulation study was started",
        "",
        "Design",
        "Baseline PSmelding and independent-clustering delta=2 Melding each used",
        "three starting areal variances: 0.001, 0.1 and 0.4; point start was 0.1",
        "All phases retain the same two saved datasets, 1676-vertex production mesh,",
        "separate observation precisions with prior B (pc.prec c(0.316,0.5)),",
        "theta/x initial vectors, seeds, 1:1 threads and 1000 latent-score draws",
        "These are controlled warm-start diagnostics, not independent cold starts",
        "",
        "One additional control changed per phase:",
        "  Original: outer tolerance 0.005, inner tolerance 0.001, internal.opt TRUE",
        "  Tight outer: outer tolerance 0.0001",
        "  Fixed backend: also internal.opt FALSE",
        "  Tight inner: also tolerance.step 0.00001",
        "",
        "Maximum differences across the three starts within each phase/model:",
        "Phase / dataset                    mean/SD   endpoint/SD   prediction RMS   log-mode spread   pass",
    ]
    for row in stability:
        label = row['phase'] + '/' + ('baseline' if row['dataset']=='baseline' else 'delta2')
        lines.append(
            f"{label:34} {float(row['maximum_parameter_mean_shift_pooled_sd']):8.5f} "
            f"{float(row['maximum_interval_endpoint_shift_pooled_sd']):12.5f} "
            f"{float(row['max_prediction_mean_RMS']):16.6f} "
            f"{float(row['log_mode_spread']):17.6f} {row['all_components_pass']}"
        )
    lines.extend([
        "",
        "Interpretation",
        "Tightening the outer tolerance clears the tested delta=2 Melding group",
        "under the recorded screens. Neither disabling internal online optimization",
        "nor tightening the inner tolerance clears baseline PSmelding. The final",
        "inner-tolerance phase reproduces the preceding deterministic summaries",
        "Its maximum mean shift is 0.20025 pooled SD, endpoint shift 0.64522 pooled",
        "SD and prediction RMS 0.010871. Theta's upper endpoint is most sensitive",
        "All 18 new fits and scores completed with mode status 0 and no negative",
        "Hessian eigenvalues; this does not establish stable posterior uncertainty",
        "No new production numerical settings were adopted or frozen",
        "A focused Hessian/finite-difference and integration investigation remains",
        "appropriate before another readiness assessment; no further fits were run",
        "",
        "The existing screens were retained: mean/median/interval differences <=0.1",
        "pooled posterior SD, latent prediction mean RMS <=0.01, log-mode spread",
        "<=0.01, plus parameter validity, optimizer, Hessian, integration and boundary",
        "screens. These are practical review thresholds, not guarantees of accuracy",
        "Means and interval endpoints come from retained posterior marginals and",
        "predictions from INLA summaries; failures cannot be fixed by more CRPS draws",
        "Integration excess is read from logs rounded to three decimals; passing at",
        "0.001 is a pass of the recorded screen rather than full-precision proof",
        "Boundary probabilities must be read alongside each phase's support audit:",
        "a zero probability beyond tabulated marginal support is not exact tail mass",
        "",
        "Scores and preservation",
        "all_phase_attempts.csv records every attempt, runtime, MSE and latent CRPS",
        "score_ranges.csv summarizes the observed MSE/CRPS ranges within each phase",
        "CRPS_Z_joint uses 1000 posterior draws; block MCSE is an approximation",
        "phase_stability.csv records every screening component and score spread",
        "Each phase directory retains raw/full fits, data, posterior draws, pointwise",
        "predictions, checkpoints, scripts, controls audits and SHA256 manifests",
        "Scripts and RDS provenance keep their original /tmp execution paths; each",
        "phase's path_mapping.csv explains access through the durable archive",
        "",
        "Laptop 2",
        "macOS with Apple silicon can target the same R architecture; its actual",
        "R/packages/native libraries and numerical results remain unverified",
        "laptop_2_apple_silicon.txt gives a no-fit environment check and limitations",
        "The current main-code laptop reference still uses default tolerance 0.005",
        "After numerical settings are accepted, a new frozen configuration and fresh",
        "laptop-1 reference are required before a fresh independent laptop-2 check",
        "Never use copied completed fits as laptop 2's output directory",
        "",
        "Publication",
        "The existing production_v2_20261022 frozen configuration stays unchanged",
        "Only this follow-up's compact reports are published to its analysis folder",
        "Full fitted objects are retained locally; the reports ZIP does not contain them",
    ])
    (ROOT / "findings.txt").write_text("\n".join(lines) + "\n")
