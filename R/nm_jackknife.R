#  ---------------------------------------------------------------------------
#  nm_jackknife.R  --  Leave-one-out (subject) jackknife for a NONMEM run
#  Refactor of "SUPPLEMENTARY SCRIPT S1.07" (run_JackKnife_Supplement_v7.R)
#  into a single exported function, on the same pattern as nm_bootstrap().
#
#  Exported :  nm_jackknife()
#  S3       :  print.nm_jackknife(), plot.nm_jackknife()
#  Helpers  :  .jk_*  (dot/jk-prefixed so that sourcing every script in the
#              user's scripts_dir alongside nm_bootstrap.R causes no name clash;
#              logic is kept VERBATIM from the verified v7 script).
#
#  Method (unchanged from v7): remove one subject (ID) at a time, re-run NONMEM,
#  and measure that subject's INFLUENCE:
#    * primary  : |dOFV| = |OFV(-i) - OFV(ref)|      (likelihood influence);
#    * secondary: relative THETA shifts (th(-i)-th_ref)/th_ref, their L2 norm
#                 ||dtheta||, the max |shift|, and an SE-standardised L2 norm
#                 ||(th(-i)-th_ref)/SE_ref|| when the reference $COV is usable.
#  Subjects are ranked by |dOFV|; two CSVs (all runs + successful summary) and
#  three figures are written to output_dir, exactly as in v7.
#  ---------------------------------------------------------------------------

library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(tibble)

# Signed integer/decimal with optional scientific exponent
.JK_NUM_PAT <- "-?[0-9]*\\.?[0-9]+(?:[Ee][+-]?[0-9]+)?"

# ===========================================================================
#  Helpers (verbatim logic from v7, renamed .jk_*)
# ===========================================================================

.jk_ensure_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  x
}

.jk_first_num <- function(x) {
  m <- regexpr(.JK_NUM_PAT, x, perl = TRUE)
  if (m < 0) return(NA_real_)
  as.numeric(regmatches(x, m))
}

.jk_extract_data_path <- function(ctl) {
  lines <- readLines(ctl)
  i <- grep("^\\s*\\$DATA", lines, ignore.case = TRUE)[1]
  if (is.na(i)) stop("No $DATA record found in control file: ", ctl)
  l <- sub("^\\s*\\$DATA\\s+", "", lines[i], ignore.case = TRUE)
  l <- strsplit(l, ";")[[1]][1]
  l <- trimws(l)
  if (grepl('^["\']', l)) {
    q <- substr(l, 1, 1)
    f <- sub(paste0("^", q, "([^", q, "]*)", q, ".*$"), "\\1", l)
  } else {
    f <- strsplit(l, "\\s+")[[1]][1]
  }
  if (!grepl("^([A-Za-z]:|/|\\\\)", f)) f <- file.path(dirname(ctl), f)
  normalizePath(f, winslash = "/", mustWork = FALSE)
}

.jk_extract_input_names <- function(ctl) {
  lines <- readLines(ctl)
  i <- grep("^\\s*\\$INPUT", lines, ignore.case = TRUE)[1]
  if (is.na(i)) return(character(0))
  block <- lines[i]; k <- i + 1
  while (k <= length(lines) && !grepl("^\\s*\\$", lines[k])) {
    block <- c(block, lines[k]); k <- k + 1
  }
  txt  <- paste(block, collapse = " ")
  txt  <- sub("^\\s*\\$INPUT\\s*", "", txt, ignore.case = TRUE)
  txt  <- strsplit(txt, ";")[[1]][1]
  toks <- strsplit(trimws(txt), "\\s+")[[1]]
  toks <- toks[nchar(toks) > 0]
  if (length(toks) == 0) return(character(0))
  nm <- vapply(toks, function(t) {
    parts <- strsplit(t, "=")[[1]]
    lab   <- parts[!toupper(parts) %in% c("DROP", "SKIP")]
    if (length(lab)) lab[1] else "DROP"
  }, character(1))
  make.unique(nm, sep = "_")
}

.jk_patch_data_line <- function(line, new_file) {
  m <- regmatches(line, regexec('^(\\s*\\$DATA\\s+)("[^"]*"|\'[^\']*\'|\\S+)(.*)$',
                                line, ignore.case = TRUE))[[1]]
  if (length(m) == 0) return(paste0("$DATA ", new_file))
  paste0(m[2], new_file, m[4])
}

.jk_strip_record <- function(ctl_lines, record) {
  pat  <- paste0("^\\s*\\$", record)
  drop <- rep(FALSE, length(ctl_lines))
  i <- 1
  while (i <= length(ctl_lines)) {
    if (grepl(pat, ctl_lines[i], ignore.case = TRUE)) {
      drop[i] <- TRUE; j <- i + 1
      while (j <= length(ctl_lines) && !grepl("^\\s*\\$", ctl_lines[j])) {
        drop[j] <- TRUE; j <- j + 1
      }
      i <- j
    } else i <- i + 1
  }
  ctl_lines[!drop]
}

.jk_extract_ofv <- function(lst) {
  lines <- readLines(lst, warn = FALSE)
  idx <- grep("#OBJV", lines)
  if (length(idx) == 0) return(NA_real_)
  .jk_first_num(lines[idx[length(idx)]])
}

.jk_parse_theta_vector <- function(lst_lines, se = FALSE) {
  if (se) {
    h0 <- grep("STANDARD ERROR OF ESTIMATE", lst_lines)
    if (length(h0) == 0) return(numeric(0))
    th <- grep("THETA - VECTOR OF FIXED EFFECTS PARAMETERS", lst_lines)
    th <- th[th > h0[1]]
  } else {
    th <- grep("THETA - VECTOR OF FIXED EFFECTS PARAMETERS", lst_lines)
  }
  if (length(th) == 0) return(numeric(0))
  start <- th[1] + 1
  omega <- grep("OMEGA", lst_lines); omega <- omega[omega > start]
  stop_line <- if (length(omega)) omega[1] - 1 else length(lst_lines)
  vals <- numeric(0)
  for (ln in lst_lines[start:stop_line]) {
    if (!nzchar(trimws(ln))) next
    if (grepl("\\bTH\\s*[0-9]", ln)) next
    for (t in strsplit(trimws(ln), "\\s+")[[1]]) {
      vals <- c(vals,
                if (grepl(paste0("^", .JK_NUM_PAT, "$"), t, perl = TRUE)) as.numeric(t)
                else NA_real_)
    }
  }
  vals
}

.jk_theta_metrics <- function(th, ref_theta, ref_se = NULL) {
  k   <- length(ref_theta)
  out <- list(rel = rep(NA_real_, k), norm = NA_real_, max = NA_real_,
              z_norm = NA_real_)
  if (k == 0 || length(th) != k) return(out)
  rel <- ifelse(ref_theta != 0, (th - ref_theta) / ref_theta, NA_real_)
  out$rel <- rel
  fin <- is.finite(rel)
  if (any(fin)) {
    out$norm <- sqrt(sum(rel[fin]^2))
    out$max  <- max(abs(rel[fin]))
  }
  if (!is.null(ref_se) && length(ref_se) == k) {
    z  <- ifelse(is.finite(ref_se) & ref_se > 0, (th - ref_theta) / ref_se, NA_real_)
    zf <- is.finite(z)
    if (any(zf)) out$z_norm <- sqrt(sum(z[zf]^2))
  }
  out
}

.jk_extract_status <- function(lst) {
  txt           <- paste(readLines(lst, warn = FALSE), collapse = "\n")
  min_ok        <- grepl("MINIMIZATION SUCCESSFUL", txt)
  cov_requested <- grepl("COVARIANCE STEP OMITTED:\\s*NO", txt)
  cov_matrix_ok <- grepl("COVARIANCE MATRIX OF ESTIMATE", txt)
  cov_singular  <- grepl("MATRIX ALGORITHMICALLY SINGULAR|R MATRIX SINGULAR|S MATRIX SINGULAR", txt)
  cov_ok        <- cov_requested & cov_matrix_ok & !cov_singular
  tibble(min_ok = min_ok, cov_ok = cov_ok, run_ok = min_ok)
}

# ===========================================================================
#  Exported function
# ===========================================================================

#' Leave-one-out (subject) jackknife of a NONMEM model
#'
#' @param ctl            Path to the reference control stream (.ctl).
#' @param nmfe_path      Path to the nmfe*.bat launcher.
#' @param lst            Reference output (.lst). Default: same basename as ctl, .lst.
#' @param reference_data Data file. NULL = auto-detected from $DATA in ctl.
#' @param id_col         Subject ID column name (default "ID").
#' @param data_has_header TRUE if the data file has an uncommented header row.
#' @param comment_char   Data-file comment char (default "#"); "" = none.
#' @param default_colnames Fallback column names if neither a commented header
#'                       nor $INPUT can supply them (default NULL).
#' @param output_dir     Working directory (default "<ctl_dir>/<base>_jackknife").
#' @param overwrite_existing_runs FALSE = resume (skip runs already minimised).
#' @param strip_tables   Remove $TABLE from the jackknife controls (default TRUE).
#' @param strip_cov      Remove $COVARIANCE from the N replicate runs (default FALSE,
#'                       faithful to v7; TRUE is faster and only affects per-run cov_ok,
#'                       NOT the |dOFV|/dtheta influence metrics).
#' @param dtheta_threshold Relative shift flagged as "large" (default 0.01).
#' @param run            FALSE = only build datasets + controls (no NONMEM).
#' @param verbose        Progress messages (default TRUE).
#'
#' @return An object of class "nm_jackknife" (see print()/plot()).
#' @export
nm_jackknife <- function(ctl,
                         nmfe_path               = NULL,
                         lst                     = NULL,
                         reference_data          = NULL,
                         id_col                  = "ID",
                         data_has_header         = FALSE,
                         comment_char            = "#",
                         default_colnames        = NULL,
                         output_dir              = NULL,
                         overwrite_existing_runs = FALSE,
                         strip_tables            = TRUE,
                         strip_cov               = FALSE,
                         dtheta_threshold        = 0.01,
                         run                     = TRUE,
                         verbose                 = TRUE) {

  say <- function(...) if (isTRUE(verbose)) message(...)

  stopifnot(file.exists(ctl))
  reference_ctl <- normalizePath(ctl, winslash = "/", mustWork = TRUE)
  ctl_dir <- dirname(reference_ctl)
  base    <- tools::file_path_sans_ext(basename(reference_ctl))
  if (is.null(lst)) lst <- file.path(ctl_dir, paste0(base, ".lst"))
  reference_lst <- normalizePath(lst, winslash = "/", mustWork = FALSE)
  if (isTRUE(run) && !file.exists(reference_lst))
    stop("Reference .lst not found: ", reference_lst)
  if (is.null(output_dir)) output_dir <- file.path(ctl_dir, paste0(base, "_jackknife"))

  # --- Directories ---------------------------------------------------------
  .jk_ensure_dir(output_dir)
  jk_data_dir <- .jk_ensure_dir(file.path(output_dir, "jk_data"))
  jk_ctl_dir  <- .jk_ensure_dir(file.path(output_dir, "jk_ctl"))
  jk_run_dir  <- .jk_ensure_dir(file.path(output_dir, "jk_run"))
  fig_dir     <- .jk_ensure_dir(file.path(output_dir, "figures"))

  # --- Read data -----------------------------------------------------------
  if (is.null(reference_data)) reference_data <- .jk_extract_data_path(reference_ctl)
  input_names <- .jk_extract_input_names(reference_ctl)

  commented_names <- NULL
  if (nchar(comment_char) > 0 && !data_has_header) {
    raw_lines  <- readLines(reference_data, n = 5)
    first_line <- raw_lines[nchar(trimws(raw_lines)) > 0][1]
    if (!is.na(first_line) && grepl(paste0("^", comment_char), trimws(first_line))) {
      header_line     <- sub(paste0("^", comment_char), "", trimws(first_line))
      commented_names <- trimws(strsplit(header_line, ",")[[1]])
    }
  }

  dat <- readr::read_csv(reference_data,
                         col_names      = data_has_header,
                         comment        = if (nchar(comment_char) > 0) comment_char else "",
                         show_col_types = FALSE)

  if (!data_has_header) {
    if (!is.null(commented_names)) {
      chosen <- commented_names; src <- "commented header"
    } else if (length(input_names) > 0) {
      chosen <- input_names; src <- "$INPUT"
    } else if (!is.null(default_colnames)) {
      chosen <- default_colnames; src <- "default_colnames"
    } else {
      stop("No column names available (no commented header, no $INPUT, ",
           "no default_colnames).")
    }
    if (length(chosen) != ncol(dat)) {
      stop("Column-name / column-count mismatch: ", length(chosen),
           " names from ", src, " but the data has ", ncol(dat), " columns.\n",
           "Names: ", paste(chosen, collapse = ", "))
    }
    colnames(dat) <- chosen
    say("Column names assigned from ", src, ": ", paste(chosen, collapse = ", "))
  }

  if (!(id_col %in% colnames(dat)))
    stop("ID column '", id_col, "' not found. Available columns: ",
         paste(colnames(dat), collapse = ", "))

  dat <- dat %>% mutate(.ID = .data[[id_col]])
  ids <- sort(unique(dat$.ID))
  N   <- length(ids)
  say(nrow(dat), " rows, ", N, " subjects: ", paste(ids, collapse = ", "))

  # --- Generate leave-one-out datasets + controls --------------------------
  base_ctl <- readLines(reference_ctl)
  if (strip_tables) base_ctl <- .jk_strip_record(base_ctl, "TABLE")
  # strip_cov : le SE de reference (distance z) provient du .lst de reference deja
  # tourne ; les runs jackknife n'ont pas besoin de $COV. Le retirer accelere sans
  # affecter |dOFV|/Δθ (seul cov_ok par run devient FALSE). Defaut FALSE = fidele au v7.
  if (isTRUE(strip_cov)) base_ctl <- .jk_strip_record(base_ctl, "COV")
  if (is.na(grep("^\\s*\\$DATA", base_ctl, ignore.case = TRUE)[1]))
    stop("No $DATA record found in the reference control; cannot patch it.")

  for (id in ids) {
    dat_jk <- dat %>% filter(.ID != id)
    name   <- paste0("jk_", id)
    utils::write.table(dat_jk %>% select(-.ID),
                       file = file.path(jk_data_dir, paste0(name, ".csv")),
                       sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    new_ctl      <- base_ctl
    idx          <- grep("^\\s*\\$DATA", new_ctl, ignore.case = TRUE)[1]
    new_ctl[idx] <- .jk_patch_data_line(new_ctl[idx], paste0(name, ".csv"))
    writeLines(new_ctl, file.path(jk_ctl_dir, paste0(name, ".ctl")))
  }
  say(N, " jackknife datasets and control files generated.")

  if (!isTRUE(run)) {
    say("run = FALSE: datasets + controls generated, no NONMEM launched.")
    return(structure(list(call = match.call(), reference_ctl = reference_ctl,
                          reference_lst = reference_lst, reference_data = reference_data,
                          output_dir = output_dir, n_subjects = N, ids = ids,
                          ref_ofv = NA_real_, ref_theta = numeric(0), ref_se = numeric(0),
                          k_theta = 0L, all_runs = NULL, summary = NULL,
                          dtheta_threshold = dtheta_threshold, run = FALSE),
                     class = "nm_jackknife"))
  }

  if (is.null(nmfe_path) || !file.exists(nmfe_path))
    stop("A valid `nmfe_path` is required when run = TRUE.")
  nmfe_win <- normalizePath(nmfe_path, winslash = "\\", mustWork = FALSE)

  # --- Reference run: OFV, THETA, SE, status -------------------------------
  ref_lines  <- readLines(reference_lst, warn = FALSE)
  ref_status <- .jk_extract_status(reference_lst)
  ref_ofv    <- .jk_extract_ofv(reference_lst)
  ref_theta  <- .jk_parse_theta_vector(ref_lines, se = FALSE)
  ref_se     <- .jk_parse_theta_vector(ref_lines, se = TRUE)
  k_theta    <- length(ref_theta)
  say("Reference OFV: ", ref_ofv)
  if (is.na(ref_ofv))
    stop("Could not read the reference OFV (#OBJV) from ", reference_lst, ". Aborting.")
  if (!isTRUE(ref_status$min_ok))
    warning("Reference run did NOT minimise successfully; dOFV are relative to a non-converged reference.")
  if (k_theta == 0)
    warning("No THETA vector read from the reference .lst; dtheta metrics will be NA.")
  if (length(ref_se) != k_theta || all(is.na(ref_se))) {
    say("[INFO] No usable THETA SE in the reference run: SE-standardised distance will be NA.")
    ref_se <- rep(NA_real_, k_theta)
  }

  # --- Run NONMEM (sequential; resume-aware) -------------------------------
  old_wd <- setwd(jk_run_dir); on.exit(setwd(old_wd), add = TRUE)
  for (id in ids) {
    name     <- paste0("jk_", id)
    lst_file <- paste0(name, ".lst")
    if (!overwrite_existing_runs && file.exists(lst_file) &&
        isTRUE(.jk_extract_status(lst_file)$min_ok)) {
      say("Skipping (already completed successfully): ", name)
      next
    }
    file.copy(file.path(jk_ctl_dir,  paste0(name, ".ctl")), ".", overwrite = TRUE)
    file.copy(file.path(jk_data_dir, paste0(name, ".csv")), ".", overwrite = TRUE)
    say("Running: ", name)
    status <- system2("cmd.exe",
                      args = c("/c", shQuote(nmfe_win, type = "cmd"),
                               paste0(name, ".ctl"), lst_file))
    if (!identical(as.integer(status), 0L))
      warning("NONMEM returned a non-zero exit status (", status, ") for ", name, ".")
  }

  # --- Extract results -----------------------------------------------------
  lst_files <- list.files(jk_run_dir, pattern = "^jk_[0-9]+\\.lst$", full.names = TRUE)
  if (length(lst_files) == 0)
    stop("No valid .lst files found in ", jk_run_dir, ".")
  say(length(lst_files), " valid .lst files found.")

  res_all <- purrr::map_dfr(lst_files, function(f) {
    id_val <- sub("^jk_([0-9]+)\\.lst$", "\\1", basename(f))
    th <- .jk_parse_theta_vector(readLines(f, warn = FALSE), se = FALSE)
    tm <- .jk_theta_metrics(th, ref_theta, ref_se)
    row <- tibble(file = basename(f), id_removed = id_val,
                  ofv = .jk_extract_ofv(f),
                  dtheta_norm = tm$norm, dtheta_max = tm$max, dtheta_znorm = tm$z_norm) %>%
      bind_cols(.jk_extract_status(f))
    if (k_theta > 0) {
      rel_cols <- as_tibble(setNames(as.list(tm$rel),
                                     paste0("dtheta_rel_", seq_len(k_theta))))
      row <- bind_cols(row, rel_cols)
    }
    row
  })

  # --- Build summary -------------------------------------------------------
  res_all <- res_all %>%
    mutate(ref_ofv = ref_ofv, dOFV = ofv - ref_ofv, abs_dOFV = abs(dOFV),
           large_dtheta = !is.na(dtheta_max) & dtheta_max > dtheta_threshold)
  readr::write_csv(res_all, file.path(output_dir, "jackknife_all_runs.csv"))

  n_failed <- sum(!res_all$run_ok); n_no_ofv <- sum(is.na(res_all$dOFV))
  if (n_failed > 0 || n_no_ofv > 0)
    say("[WARNING] Excluded from ranking: ", n_failed, " non-minimised and ",
        n_no_ofv, " missing-OFV run(s).")

  res <- res_all %>% filter(run_ok, is.finite(dOFV)) %>%
    arrange(desc(abs_dOFV)) %>% mutate(rank = row_number())
  if (nrow(res) == 0)
    stop("No usable run remained after filtering (all failed or missing OFV).")

  res <- res %>% arrange(desc(dtheta_norm)) %>% mutate(rank_dtheta = row_number()) %>%
    arrange(rank)
  if (any(is.finite(res$dtheta_znorm)))
    res <- res %>% arrange(desc(dtheta_znorm)) %>% mutate(rank_dtheta_z = row_number()) %>%
      arrange(rank)

  readr::write_csv(res, file.path(output_dir, "jackknife_summary.csv"))
  say("Summary written to: ", file.path(output_dir, "jackknife_summary.csv"))

  # --- Figures (written to disk, as in v7; also reproducible via plot()) ----
  obj <- structure(list(call = match.call(), reference_ctl = reference_ctl,
                        reference_lst = reference_lst, reference_data = reference_data,
                        output_dir = output_dir, n_subjects = N, ids = ids,
                        ref_ofv = ref_ofv, ref_theta = ref_theta, ref_se = ref_se,
                        k_theta = k_theta, all_runs = res_all, summary = res,
                        dtheta_threshold = dtheta_threshold, run = TRUE),
                   class = "nm_jackknife")
  try({
    ggsave(file.path(fig_dir, "abs_dOFV_bar.png"),  plot(obj, "bar"),     width = 6, height = 5)
    ggsave(file.path(fig_dir, "rank_dOFV.png"),     plot(obj, "rank"),    width = 7, height = 4)
    if (obj$k_theta > 0 && any(is.finite(res$dtheta_norm)))
      ggsave(file.path(fig_dir, "dOFV_vs_dtheta.png"), plot(obj, "scatter"), width = 6, height = 5)
  }, silent = TRUE)
  say("Figures saved to: ", fig_dir)
  say("Jackknife analysis complete.")
  obj
}

# ===========================================================================
#  S3 methods
# ===========================================================================

#' @export
print.nm_jackknife <- function(x, digits = 4, ...) {
  cat("<nm_jackknife>\n")
  cat(sprintf("  reference : %s\n", basename(x$reference_ctl)))
  cat(sprintf("  run dir   : %s\n", x$output_dir))
  cat(sprintf("  subjects  : %d\n", x$n_subjects))
  if (!isTRUE(x$run) || is.null(x$summary)) { cat("  (dry run: datasets + controls only)\n"); return(invisible(x)) }
  cat(sprintf("  ref OFV   : %s\n", signif(x$ref_ofv, digits)))
  cat(sprintf("  runs used : %d (of %d)\n", nrow(x$summary), nrow(x$all_runs)))
  top <- utils::head(x$summary[, c("rank", "id_removed", "dOFV", "abs_dOFV",
                                   "dtheta_norm", "dtheta_max")], 10)
  num <- vapply(top, is.numeric, logical(1)); top[num] <- lapply(top[num], round, digits)
  cat("  most influential subjects (by |dOFV|):\n")
  print(top, row.names = FALSE)
  invisible(x)
}

#' @export
plot.nm_jackknife <- function(x, type = c("bar", "rank", "scatter"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot.nm_jackknife().")
  res <- x$summary
  if (is.null(res) || !nrow(res)) stop("No summary to plot (dry run?).")

  if (type == "bar") {
    return(ggplot2::ggplot(res, ggplot2::aes(x = stats::reorder(factor(.data$id_removed), .data$abs_dOFV),
                                             y = .data$abs_dOFV)) +
             ggplot2::geom_col(fill = "steelblue") + ggplot2::coord_flip() +
             ggplot2::theme_bw() +
             ggplot2::labs(x = "Subject removed", y = "|\u0394OFV|",
                           title = "Jackknife \u2014 influence by subject"))
  }
  if (type == "rank") {
    return(ggplot2::ggplot(res, ggplot2::aes(x = .data$rank, y = .data$abs_dOFV,
                                             label = .data$id_removed)) +
             ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
             ggplot2::geom_text(vjust = -0.8, size = 3) + ggplot2::theme_bw() +
             ggplot2::labs(x = "Rank", y = "|\u0394OFV|",
                           title = "Jackknife \u2014 influence ranking"))
  }
  if (!(x$k_theta > 0 && any(is.finite(res$dtheta_norm))))
    stop("No THETA-shift data to plot (scatter needs a reference THETA vector).")
  ggplot2::ggplot(res, ggplot2::aes(x = .data$dtheta_norm, y = .data$abs_dOFV,
                                    label = .data$id_removed)) +
    ggplot2::geom_point(size = 2) + ggplot2::geom_text(vjust = -0.8, size = 3) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "\u2016\u0394\u03b8\u2016 (relative parameter shift)", y = "|\u0394OFV|",
                  title = "Jackknife \u2014 likelihood influence vs parameter sensitivity")
}
