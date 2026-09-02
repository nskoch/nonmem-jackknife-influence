nonmem-jackknife-influence

A single, self-contained R function — nm_jackknife() — that performs a leave-one-subject-out jackknife on a fitted NONMEM model to quantify how much each individual influences the fit.

No PsN, no xpose — base R plus the tidyverse.

This script is also bundled inside NMtree (a Shiny application to drive NONMEM end-to-end). It is published here on its own so it can be found, used, and cited independently. See [NMtree](https://github.com/<< your-account >>/nmtree).

What it does

For a reference run, nm_jackknife() removes one subject (ID) at a time, re-estimates the model in NONMEM, and measures that subject's influence:

Primary metric — |dOFV| = |OFV(-i) - OFV(ref)| (likelihood influence).
Secondary metrics — relative THETA shifts (theta(-i) - theta_ref) / theta_ref, their L2 norm ||dtheta||, the maximum |shift|, and an SE-standardised L2 norm ||(theta(-i) - theta_ref) / SE_ref|| when the reference $COVARIANCE is usable.

Subjects are ranked by |dOFV|. The function writes two CSVs (jackknife_all_runs.csv, jackknife_summary.csv) and three figures to the output directory, and returns an nm_jackknife object with print() and plot() methods.

Requirements
R (>= 4.1) with: dplyr, readr, purrr, ggplot2, tibble.
r
  install.packages(c("dplyr", "readr", "purrr", "ggplot2", "tibble"))
NONMEM installed, with its nmfe*.bat launcher — required for the re-estimation step. The re-run step is driven through cmd.exe, i.e. it is Windows-oriented.
You can run with run = FALSE to only build the N datasets and control streams (no NONMEM, cross-platform), then estimate them yourself.
Usage
r
source("nm_jackknife.R")

jk <- nm_jackknife(
  ctl       = "path/to/run.ctl",         # reference control stream
  nmfe_path = "C:/nm74/util/nmfe74.bat"  # NONMEM launcher
  # lst, reference_data, id_col, output_dir, ... are auto-detected / optional
)

print(jk)           # ranked influence table
plot(jk, "bar")     # |dOFV| per subject
plot(jk, "rank")    # ranked |dOFV|
plot(jk, "scatter") # |dOFV| vs ||dtheta||
Key arguments
Argument	Default	Role
ctl	—	Reference control stream (.ctl). Required.
nmfe_path	NULL	Path to nmfe*.bat. Needed when run = TRUE.
lst	<ctl>.lst	Reference output; defaults to the ctl basename.
reference_data	NULL	Data file; auto-detected from $DATA if NULL.
id_col	"ID"	Subject identifier column.
data_has_header	FALSE	TRUE if the data file has an uncommented header.
comment_char	"#"	Data-file comment character ("" = none).
output_dir	<ctl_dir>/<base>_jackknife	Where datasets, controls, runs, CSVs and figures go.
overwrite_existing_runs	FALSE	FALSE = resume (skip already-minimised runs).
strip_tables	TRUE	Drop $TABLE from the jackknife controls.
strip_cov	FALSE	Drop $COVARIANCE from replicate runs (faster; does not affect the `
dtheta_threshold	0.01	Relative shift flagged as "large".
run	TRUE	FALSE = build datasets + controls only, no NONMEM.
verbose	TRUE	Progress messages.
Outputs

Written to output_dir:

jackknife_all_runs.csv — every replicate (including non-converged).
jackknife_summary.csv — successful runs, ranked by |dOFV|.
abs_dOFV_bar.png, rank_dOFV.png, dOFV_vs_dtheta.png.

The returned object also carries the reference OFV/THETA/SE, the per-subject tables ($all_runs, $summary), and the settings used.

Method note

|dOFV| captures how much each subject moves the objective function; the THETA shift metrics capture how much each subject moves the parameter estimates. A subject high on both is genuinely influential. Interpretation is descriptive — this is an influence diagnostic, not a formal outlier test.

License
<!-- Add a LICENSE file (Add file > Create new file > name it "LICENSE" > pick a template). MIT is a common, permissive choice; match whatever you use for NMtree. -->

(To be added — e.g. MIT.)

Citation

If you use this script, please cite it. A CITATION.cff can be added so GitHub shows a "Cite this repository" button; a Zenodo release gives it a DOI.
