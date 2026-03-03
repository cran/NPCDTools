# NPCDTools 1.1.0

## New Functions

- `distractor.check()`: Detects implausible and improper distractors in a Q-matrix for multiple-choice items (Chiu, Köhn, & Wang, in press).
- `Q.implausible()`: Generates a Q-matrix containing implausible MC item distractors from a proper and plausible Q-matrix.
- `Q.improper()`: Generates a Q-matrix containing improper MC item distractors from a proper and plausible Q-matrix.
- `plot.GNPC()`: S3 plot method for GNPC objects, providing convergence tracking visualization for individual examinees.
- `run_gnpc_app()`: Launches an interactive Shiny application demonstrating NPC, GNPC, and G-DINA workflows.

## New Dataset

- `Q_Ozaki`: A Q-matrix for 30 multiple-choice items measuring 5 attributes from Ozaki (2015), including coded distractors.

## Improvements

- Added S3 `print` methods for `GNPC`, `NPC`, `QR`, `Q.completeness`, and `distractor.check` objects, providing structured and informative summaries.
- Added convergence tracking option in `GNPC()` via the `track.convergence` argument.
- Improved input validation across functions via internal `CheckInput()`.
- Consolidated internal pattern generation for consistency across all functions.

## Dependency Changes

- Removed dependency on `NPCD`. The `TSQE()` function now uses the package's own `QR()` function for Q-matrix refinement instead of `NPCD::Qrefine()`.
- Removed dependency on `SimDesign`.

# NPCDTools 1.0

- Initial CRAN release (2024-09-23).
- Core nonparametric classification methods: `NPC()` and `GNPC()`.
- Q-matrix tools: `Q.completeness()`, `Q.generate()`, `QR()`, `TSQE()`, `bestQperm()`.
- Evaluation metrics: `AAR()`, `PAR()`, `RR()`, `correction.rate()`, `retention.rate()`.
