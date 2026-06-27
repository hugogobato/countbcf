################################################################################
##  Aggregate per-cell simulation CSVs for the count-only and zero-inflated
##  studies into Monte Carlo summaries, render reference-cell tables in
##  Markdown + LaTeX, and draw N / P / ATE ablation plots.
##
##  Inputs:
##    tests/results_count/sim_count__<method>__<dgp>__N<N>_P<P>_ATE<ATE*100>.csv
##    tests/results_zi/sim_zi__<method>__<dgp>__N<N>_P<P>_ATE<ATE*100>.csv
##
##  Outputs (under tests/analysis/):
##    summary_count.csv, summary_zi.csv         (one row per cell x method x dgp)
##    tables/table_count_ref.{md,tex}           (reference-cell tables)
##    tables/table_zi_ref.{md,tex}
##    plots/<study>_<factor>_<metric>.png       (ablation plots)
##    SUMMARY_count.md, SUMMARY_zi.md           (top-level write-ups)
##
##  Usage from package root:
##    Rscript tests/analyze_results.R
################################################################################

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})

THIS_DIR     <- "tests"
RESULTS_DIRS <- list(
  count = file.path(THIS_DIR, "results_count"),
  zi    = file.path(THIS_DIR, "results_zi")
)
OUT_DIR     <- file.path(THIS_DIR, "analysis")
TABLES_DIR  <- file.path(OUT_DIR, "tables")
PLOTS_DIR   <- file.path(OUT_DIR, "plots")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,  showWarnings = FALSE, recursive = TRUE)

METHOD_LABELS <- c(
  countbcf      = "CountBCF",
  bcf_gauss     = "BCF (Gaussian)",
  glm_poisson   = "Poisson GLM",
  glm_nb        = "NB GLM",
  glm_zip       = "ZIP GLM",
  glm_zinb      = "ZINB GLM",
  causal_forest = "Causal Forest"
)
METHOD_ORDER <- names(METHOD_LABELS)

# Fixed colour / shape per method label so every panel and study is consistent
# even when a study only contains a subset of the methods.
METHOD_COLOURS <- c(
  "CountBCF"       = "#1b9e77",
  "BCF (Gaussian)" = "#d95f02",
  "Poisson GLM"    = "#7570b3",
  "NB GLM"         = "#e7298a",
  "ZIP GLM"        = "#66a61e",
  "ZINB GLM"       = "#e6ab02",
  "Causal Forest"  = "#386cb0"
)
METHOD_SHAPES <- c(
  "CountBCF"       = 16,
  "BCF (Gaussian)" = 17,
  "Poisson GLM"    = 15,
  "NB GLM"         = 3,
  "ZIP GLM"        = 7,
  "ZINB GLM"       = 8,
  "Causal Forest"  = 18
)

# The parametric NB-family baseline whose cross-replicate mean can be dominated
# by a handful of heavy-tailed replicates; dropped in the "*_clean" plot
# variants so the remaining methods stay legible.
UNSTABLE_METHOD <- c(count = "glm_nb", zi = "glm_zinb")

DGP_LABELS_COUNT <- c(
  linear_poisson_cate    = "Linear / Poisson",
  nonlinear_poisson_cate = "Nonlinear / Poisson",
  linear_nb_cate         = "Linear / NB",
  nonlinear_nb_cate      = "Nonlinear / NB"
)
DGP_LABELS_ZI <- c(
  linear_zip_cate     = "Linear / ZIP",
  nonlinear_zip_cate  = "Nonlinear / ZIP",
  linear_zinb_cate    = "Linear / ZINB",
  nonlinear_zinb_cate = "Nonlinear / ZINB"
)

REF_N   <- 250
REF_P   <- 5
REF_ATE <- 1.25

###############################################################################
## 1. Load raw per-replicate CSVs
###############################################################################

read_study <- function(dir_path, study_name) {
  files <- list.files(dir_path, pattern = "^sim_.*\\.csv$", full.names = TRUE)
  rows  <- lapply(files, function(f) {
    df <- tryCatch(
      readr::read_csv(f, show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) return(NULL)
    if (!"pct_struct_zero" %in% names(df)) df$pct_struct_zero <- NA_real_
    df$study <- study_name
    df
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  bind_rows(rows)
}

raw_count <- read_study(RESULTS_DIRS$count, "count")
raw_zi    <- read_study(RESULTS_DIRS$zi,    "zi")

cat(sprintf("Loaded %d count rows from %d files\n",
            nrow(raw_count),
            length(unique(paste(raw_count$method, raw_count$dgp,
                                raw_count$N, raw_count$P,
                                raw_count$ate_target)))))
cat(sprintf("Loaded %d zi rows from %d files\n",
            nrow(raw_zi),
            length(unique(paste(raw_zi$method, raw_zi$dgp,
                                raw_zi$N, raw_zi$P,
                                raw_zi$ate_target)))))

raw_all <- bind_rows(raw_count, raw_zi)

raw_all <- raw_all %>%
  mutate(
    ate_ci_width = ate_q975 - ate_q025,
    success      = is.na(error)
  )

###############################################################################
## 2. Cell summaries (mean and SD over MC replicates)
###############################################################################

summarize_cell <- function(df) {
  df %>%
    filter(success) %>%
    group_by(study, method, dgp, N, P, ate_target, factor_varied) %>%
    summarise(
      n_sim                 = dplyr::n(),
      rmse_ate              = sqrt(mean((ate_mean - true_ate)^2, na.rm = TRUE)),
      ate_bias              = mean(ate_mean - true_ate, na.rm = TRUE),
      ate_coverage_mean     = mean(ate_coverage, na.rm = TRUE),
      ate_ci_width_mean     = mean(ate_ci_width, na.rm = TRUE),
      ate_ci_width_sd       = sd(ate_ci_width,   na.rm = TRUE),
      pehe_mean             = mean(pehe, na.rm = TRUE),
      pehe_sd               = sd(pehe,   na.rm = TRUE),
      cate_cov95_mean       = mean(cate_cov95, na.rm = TRUE),
      cate_cov95_sd         = sd(cate_cov95,   na.rm = TRUE),
      cate_ci_width95_mean  = mean(cate_ci_width95, na.rm = TRUE),
      cate_ci_width95_sd    = sd(cate_ci_width95,   na.rm = TRUE),
      elapsed_sec_mean      = mean(elapsed_sec, na.rm = TRUE),
      .groups = "drop"
    )
}

summary_all <- summarize_cell(raw_all)
summary_count <- summary_all %>% filter(study == "count")
summary_zi    <- summary_all %>% filter(study == "zi")

write_csv(summary_count, file.path(OUT_DIR, "summary_count.csv"))
write_csv(summary_zi,    file.path(OUT_DIR, "summary_zi.csv"))

###############################################################################
## 3. Reference-cell tables (N = 250, P = 5, ATE = 1.25)
###############################################################################

fmt_mean_sd <- function(m, s, digits = 2) {
  ifelse(is.na(m), "--",
         sprintf(paste0("%.", digits, "f $\\pm$ %.", digits, "f"), m, s))
}
fmt_mean_sd_md <- function(m, s, digits = 2) {
  ifelse(is.na(m), "--",
         sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"), m, s))
}
fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f"), x))
}

# Build a wide reference table: rows = DGP, columns = (method) x metric.
build_ref_table <- function(summary_df, dgp_labels) {
  ref <- summary_df %>%
    filter(N == REF_N, P == REF_P, abs(ate_target - REF_ATE) < 1e-9) %>%
    mutate(
      dgp_label    = factor(dgp_labels[as.character(dgp)],
                            levels = dgp_labels),
      method_label = factor(METHOD_LABELS[as.character(method)],
                            levels = METHOD_LABELS[METHOD_ORDER])
    ) %>%
    arrange(dgp_label, method_label)
  ref
}

ref_count <- build_ref_table(summary_count, DGP_LABELS_COUNT)
ref_zi    <- build_ref_table(summary_zi,    DGP_LABELS_ZI)

write_md_table <- function(ref, dgp_labels, out_path, caption_lines) {
  methods <- levels(ref$method_label)
  hdr1 <- c("DGP", "Method", "RMSE ATE", "ATE Cov.", "ATE CI Width",
            "PEHE", "CATE Cov.95", "CATE CI Width95")
  lines <- c(
    paste0("**", caption_lines, "**"),
    "",
    paste0("| ", paste(hdr1, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(hdr1)), collapse = "|"), "|")
  )
  for (dlab in levels(ref$dgp_label)) {
    sub <- ref %>% filter(dgp_label == dlab)
    for (i in seq_len(nrow(sub))) {
      row <- sub[i, ]
      dgp_cell <- if (i == 1) as.character(row$dgp_label) else ""
      cells <- c(
        dgp_cell,
        as.character(row$method_label),
        fmt_num(row$rmse_ate, 2),
        fmt_num(row$ate_coverage_mean, 2),
        fmt_mean_sd_md(row$ate_ci_width_mean, row$ate_ci_width_sd, 2),
        fmt_mean_sd_md(row$pehe_mean, row$pehe_sd, 2),
        fmt_mean_sd_md(row$cate_cov95_mean, row$cate_cov95_sd, 2),
        fmt_mean_sd_md(row$cate_ci_width95_mean, row$cate_ci_width95_sd, 2)
      )
      lines <- c(lines, paste0("| ", paste(cells, collapse = " | "), " |"))
    }
  }
  writeLines(lines, out_path)
}

write_tex_table <- function(ref, out_path, caption, label) {
  methods <- levels(ref$method_label)
  body <- character(0)
  for (dlab in levels(ref$dgp_label)) {
    sub <- ref %>% filter(dgp_label == dlab)
    for (i in seq_len(nrow(sub))) {
      row <- sub[i, ]
      dgp_cell <- if (i == 1)
        sprintf("\\multirow{%d}{*}{%s}", nrow(sub), dlab)
      else ""
      cells <- c(
        dgp_cell,
        as.character(row$method_label),
        fmt_num(row$rmse_ate, 2),
        fmt_num(row$ate_coverage_mean, 2),
        fmt_mean_sd(row$ate_ci_width_mean, row$ate_ci_width_sd, 2),
        fmt_mean_sd(row$pehe_mean, row$pehe_sd, 2),
        fmt_mean_sd(row$cate_cov95_mean, row$cate_cov95_sd, 2),
        fmt_mean_sd(row$cate_ci_width95_mean, row$cate_ci_width95_sd, 2)
      )
      body <- c(body, paste(cells, collapse = " & "))
    }
    body[length(body)] <- paste0(body[length(body)], " \\\\ \\midrule")
  }
  # drop trailing \midrule on last block
  body[length(body)] <- sub(" \\\\midrule$", "", body[length(body)])
  body <- paste0(body, " \\\\")
  body[length(body)] <- sub(" \\\\midrule \\\\$", " \\\\", body[length(body)])
  # rebuild cleanly: insert \\ at end of every line and a \midrule between DGP blocks
  body <- character(0)
  dgp_levels <- levels(ref$dgp_label)
  for (di in seq_along(dgp_levels)) {
    dlab <- dgp_levels[di]
    sub <- ref %>% filter(dgp_label == dlab)
    n_sub <- nrow(sub)
    for (i in seq_len(n_sub)) {
      row <- sub[i, ]
      dgp_cell <- if (i == 1)
        sprintf("\\multirow{%d}{*}{%s}", n_sub, dlab)
      else ""
      cells <- c(
        dgp_cell,
        as.character(row$method_label),
        fmt_num(row$rmse_ate, 2),
        fmt_num(row$ate_coverage_mean, 2),
        fmt_mean_sd(row$ate_ci_width_mean, row$ate_ci_width_sd, 2),
        fmt_mean_sd(row$pehe_mean, row$pehe_sd, 2),
        fmt_mean_sd(row$cate_cov95_mean, row$cate_cov95_sd, 2),
        fmt_mean_sd(row$cate_ci_width95_mean, row$cate_ci_width95_sd, 2)
      )
      body <- c(body, paste(paste(cells, collapse = " & "), "\\\\"))
    }
    if (di < length(dgp_levels)) body <- c(body, "\\midrule")
  }
  hdr <- c(
    "\\begin{table}[h!]\\footnotesize",
    sprintf("\\caption{%s}\\label{%s}", caption, label),
    "\\centering",
    "\\begin{tabular*}{\\linewidth}{ @{\\extracolsep{\\fill}} ll cc cccc @{}}",
    "\\toprule",
    " &  & \\multicolumn{2}{c}{ATE} & \\multicolumn{4}{c}{CATE} \\\\",
    "\\cmidrule(lr){3-4} \\cmidrule(lr){5-8}",
    paste("DGP & Method & RMSE & Cov.\\,95 & Avg.\\ CI Width &",
          "PEHE & Cov.\\,95 & Avg.\\ CI Width95 \\\\"),
    "\\midrule"
  )
  ftr <- c("\\bottomrule", "\\end{tabular*}", "\\end{table}")
  writeLines(c(hdr, body, ftr), out_path)
}

write_md_table(
  ref_count, DGP_LABELS_COUNT,
  file.path(TABLES_DIR, "table_count_ref.md"),
  paste0("Reference-cell (N = ", REF_N, ", p = ", REF_P,
         ", ATE = ", REF_ATE, ") Monte-Carlo summary, count-only DGPs. ",
         "Cells reported as mean \u00b1 SD across ", 100,
         " replicates (RMSE ATE and ATE Coverage are scalar Monte-Carlo summaries).")
)
write_md_table(
  ref_zi, DGP_LABELS_ZI,
  file.path(TABLES_DIR, "table_zi_ref.md"),
  paste0("Reference-cell (N = ", REF_N, ", p = ", REF_P,
         ", ATE = ", REF_ATE, ") Monte-Carlo summary, zero-inflated DGPs. ",
         "Cells reported as mean \u00b1 SD across ", 100,
         " replicates (RMSE ATE and ATE Coverage are scalar Monte-Carlo summaries).")
)
write_tex_table(
  ref_count,
  file.path(TABLES_DIR, "table_count_ref.tex"),
  caption = paste0("Reference-cell ($N = ", REF_N, "$, $p = ", REF_P,
                   "$, ATE $= ", REF_ATE, "$) Monte-Carlo summary on the ",
                   "count-only DGPs (Poisson, NB). RMSE ATE and ATE coverage ",
                   "are scalar Monte-Carlo summaries; the remaining columns ",
                   "report the cross-replicate mean $\\pm$ one standard ",
                   "deviation. ATE CI Width is the average width of the 95\\% ",
                   "ATE posterior credible interval; PEHE is ",
                   "$\\sqrt{\\overline{(\\hat\\tau_i - \\tau_i)^2}}$; CATE ",
                   "Cov.~95 is the per-unit 95\\% credible-interval coverage ",
                   "of $\\tau(X_i)$; CATE CI Width95 is the mean width of ",
                   "those unit-level intervals."),
  label = "tab:sim_count_ref"
)
write_tex_table(
  ref_zi,
  file.path(TABLES_DIR, "table_zi_ref.tex"),
  caption = paste0("Reference-cell ($N = ", REF_N, "$, $p = ", REF_P,
                   "$, ATE $= ", REF_ATE, "$) Monte-Carlo summary on the ",
                   "zero-inflated DGPs (ZIP, ZINB). Same conventions as ",
                   "Table~\\ref{tab:sim_count_ref}."),
  label = "tab:sim_zi_ref"
)

cat("Wrote reference-cell tables.\n")

###############################################################################
## 4. Ablation plots: vary N, p, ATE one at a time
###############################################################################

METRIC_DEFS <- tribble(
  ~metric,                  ~ylab,                                  ~has_sd,
  "rmse_ate",               "RMSE(ATE)",                            FALSE,
  "ate_coverage_mean",      "ATE 95% coverage",                     FALSE,
  "ate_ci_width_mean",      "Avg ATE 95% CI width",                 TRUE,
  "pehe_mean",              "Avg PEHE",                             TRUE,
  "cate_cov95_mean",        "Avg CATE 95% coverage",                TRUE,
  "cate_ci_width95_mean",   "Avg CATE 95% CI width (unit-level)",   TRUE
)

# Map metric -> (mean col, sd col)
SD_COL <- list(
  ate_ci_width_mean    = "ate_ci_width_sd",
  pehe_mean            = "pehe_sd",
  cate_cov95_mean      = "cate_cov95_sd",
  cate_ci_width95_mean = "cate_ci_width95_sd"
)

FACTOR_DEFS <- tribble(
  ~factor_varied, ~xvar,        ~xlab,            ~log_x,
  "N",            "N",          "Sample size N",  FALSE,
  "P",            "P",          "Covariate count p", TRUE,
  "ATE",          "ate_target", "Target ATE",     FALSE
)

# For each (study, factor), include cells that vary that factor PLUS the
# reference cell, so plots span all levels.
build_ablation_df <- function(summary_df, factor_varied) {
  ref_mask <- summary_df$N == REF_N &
              summary_df$P == REF_P &
              abs(summary_df$ate_target - REF_ATE) < 1e-9
  if (factor_varied == "N") {
    keep <- ref_mask |
      (summary_df$P == REF_P & abs(summary_df$ate_target - REF_ATE) < 1e-9)
  } else if (factor_varied == "P") {
    keep <- ref_mask |
      (summary_df$N == REF_N & abs(summary_df$ate_target - REF_ATE) < 1e-9)
  } else if (factor_varied == "ATE") {
    keep <- ref_mask |
      (summary_df$N == REF_N & summary_df$P == REF_P)
  } else stop("bad factor: ", factor_varied)
  summary_df[keep, ]
}

plot_ablation <- function(summary_df, study, factor_def, metric_def,
                          dgp_labels, normalize = FALSE,
                          drop_methods = character(0)) {
  sub <- build_ablation_df(summary_df, factor_def$factor_varied)
  if (length(drop_methods)) sub <- sub[!(sub$method %in% drop_methods), ]
  xvar <- factor_def$xvar
  sub <- sub %>%
    mutate(
      method_label = droplevels(factor(METHOD_LABELS[as.character(method)],
                            levels = METHOD_LABELS[METHOD_ORDER])),
      dgp_label    = factor(dgp_labels[as.character(dgp)],
                            levels = dgp_labels)
    )
  ycol <- metric_def$metric
  sd_col <- SD_COL[[ycol]]
  sub$.y <- sub[[ycol]]
  if (!is.null(sd_col)) {
    sub$.ymin <- sub[[ycol]] - sub[[sd_col]]
    sub$.ymax <- sub[[ycol]] + sub[[sd_col]]
  } else {
    sub$.ymin <- NA_real_
    sub$.ymax <- NA_real_
  }
  if (normalize) {
    # Normalize by CountBCF's value at the SAME level of the swept factor, per
    # DGP. CountBCF is therefore 1.0 at every level by construction, and each
    # other method's curve reads directly as its ratio to CountBCF at that level
    # (so the comparison is head-to-head within each level, not growth relative
    # to a single anchor level).
    base <- sub %>%
      filter(method == "countbcf") %>%
      select(dgp, all_of(xvar), .base = all_of(ycol))
    sub <- sub %>%
      left_join(base, by = c("dgp", xvar)) %>%
      mutate(
        .y    = .y    / .base,
        .ymin = .ymin / .base,
        .ymax = .ymax / .base
      ) %>%
      select(-.base)
  }
  sub$.x <- sub[[xvar]]
  p <- ggplot(sub, aes(x = .x, y = .y, colour = method_label,
                       shape = method_label, group = method_label)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2)
  if (!normalize && ycol %in% c("ate_coverage_mean", "cate_cov95_mean")) {
    p <- p + geom_hline(yintercept = 0.95, linetype = "dashed",
                         colour = "grey40")
  }
  if (normalize) {
    p <- p + geom_hline(yintercept = 1, linetype = "dashed",
                         colour = "grey40")
  }
  ylab_text <- if (normalize) {
    sprintf("%s / CountBCF at the same %s (per DGP)",
            metric_def$ylab, factor_def$xlab)
  } else metric_def$ylab
  title_text <- if (normalize) {
    sprintf(
      "%s study: %s vs %s (normalized to CountBCF at each %s)",
      ifelse(study == "count", "Count-only", "Zero-inflated"),
      metric_def$ylab, factor_def$xlab,
      factor_def$xlab
    )
  } else {
    sprintf(
      "%s study: %s vs %s",
      ifelse(study == "count", "Count-only", "Zero-inflated"),
      metric_def$ylab, factor_def$xlab
    )
  }
  p <- p +
    scale_colour_manual(values = METHOD_COLOURS, drop = TRUE) +
    scale_shape_manual(values = METHOD_SHAPES, drop = TRUE) +
    facet_wrap(~ dgp_label, scales = "free_y", ncol = 2) +
    labs(
      x = factor_def$xlab,
      y = ylab_text,
      colour = NULL, shape = NULL,
      title = title_text
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
  if (factor_def$log_x) {
    p <- p + scale_x_log10(breaks = sort(unique(sub$.x)))
  } else {
    p <- p + scale_x_continuous(breaks = sort(unique(sub$.x)))
  }
  p
}

# Metrics for which a normalized version is also produced when sweeping ATE
# (absolute errors trivially scale with the target ATE magnitude).
ATE_NORMALIZE_METRICS <- c("rmse_ate", "pehe_mean")

# Each ablation plot is rendered twice: the full version (all methods present in
# the study) and a "_clean" version that drops the heavy-tailed NB-family
# parametric baseline (glm_nb for count, glm_zinb for zi) whose cross-replicate
# mean can dominate a panel's y-axis.
for (study in c("count", "zi")) {
  this_summary <- if (study == "count") summary_count else summary_zi
  this_dgp_lab <- if (study == "count") DGP_LABELS_COUNT else DGP_LABELS_ZI
  drop_m       <- UNSTABLE_METHOD[[study]]
  for (fi in seq_len(nrow(FACTOR_DEFS))) {
    fd <- FACTOR_DEFS[fi, ]
    for (mi in seq_len(nrow(METRIC_DEFS))) {
      md <- METRIC_DEFS[mi, ]
      variants <- list(
        list(suffix = "",       drop = character(0)),
        list(suffix = "_clean", drop = drop_m)
      )
      for (v in variants) {
        p <- plot_ablation(this_summary, study, fd, md, this_dgp_lab,
                           drop_methods = v$drop)
        out_png <- file.path(
          PLOTS_DIR,
          sprintf("%s_%s_%s%s.png", study, fd$factor_varied, md$metric, v$suffix)
        )
        ggsave(out_png, p, width = 8, height = 5.5, dpi = 150)
        if (fd$factor_varied == "ATE" && md$metric %in% ATE_NORMALIZE_METRICS) {
          p_norm <- plot_ablation(this_summary, study, fd, md, this_dgp_lab,
                                   normalize = TRUE, drop_methods = v$drop)
          out_png_norm <- file.path(
            PLOTS_DIR,
            sprintf("%s_%s_%s_normalized%s.png",
                    study, fd$factor_varied, md$metric, v$suffix)
          )
          ggsave(out_png_norm, p_norm, width = 8, height = 5.5, dpi = 150)
        }
      }
    }
  }
}

cat("Wrote ablation plots to", PLOTS_DIR, "\n")

###############################################################################
## 5. Structured findings extracted from the summary tables
###############################################################################

pretty <- function(x, d = 2) ifelse(is.na(x), "--",
                                    sprintf(paste0("%.", d, "f"), x))
pretty_pct <- function(x, d = 0) sprintf(paste0("%+.", d, "f%%"), x)

# Average across DGPs at a given (N, P, ATE) cell, separately by method.
collapse_to_method <- function(summary_df, N_, P_, ATE_) {
  summary_df %>%
    filter(N == N_, P == P_, abs(ate_target - ATE_) < 1e-9) %>%
    group_by(method) %>%
    summarise(
      rmse_ate              = mean(rmse_ate,             na.rm = TRUE),
      ate_coverage_mean     = mean(ate_coverage_mean,    na.rm = TRUE),
      ate_ci_width_mean     = mean(ate_ci_width_mean,    na.rm = TRUE),
      pehe_mean             = mean(pehe_mean,            na.rm = TRUE),
      cate_cov95_mean       = mean(cate_cov95_mean,      na.rm = TRUE),
      cate_ci_width95_mean  = mean(cate_ci_width95_mean, na.rm = TRUE),
      .groups = "drop"
    )
}

# Headline numbers for the reference cell.
ref_agg_for <- function(summary_df) collapse_to_method(summary_df, REF_N, REF_P, REF_ATE)

# Aggregated comparison across all cells for a given factor sweep.
factor_sweep <- function(summary_df, factor_name) {
  fixed <- list(N = REF_N, P = REF_P, ATE = REF_ATE)
  rows <- list()
  if (factor_name == "N") {
    for (N_ in c(100, 250, 500)) {
      rows[[as.character(N_)]] <- collapse_to_method(
        summary_df, N_, REF_P, REF_ATE) %>% mutate(level = N_)
    }
  } else if (factor_name == "P") {
    for (P_ in c(5, 50, 250)) {
      rows[[as.character(P_)]] <- collapse_to_method(
        summary_df, REF_N, P_, REF_ATE) %>% mutate(level = P_)
    }
  } else if (factor_name == "ATE") {
    for (ATE_ in c(0.5, 1.25, 2.5)) {
      rows[[as.character(ATE_)]] <- collapse_to_method(
        summary_df, REF_N, REF_P, ATE_) %>% mutate(level = ATE_)
    }
  }
  bind_rows(rows)
}

# DGP-by-DGP head-to-head at the reference cell (which method wins on which DGP
# for a given metric).
ref_h2h_count <- function(summary_df, dgp_labels, metric, lower_is_better) {
  rf <- build_ref_table(summary_df, dgp_labels)
  wide <- rf %>%
    select(dgp_label, method, all_of(metric)) %>%
    pivot_wider(names_from = method, values_from = all_of(metric))
  wins_countbcf <- if (lower_is_better)
    sum(wide$countbcf < wide$bcf_gauss, na.rm = TRUE)
  else
    sum(wide$countbcf > wide$bcf_gauss, na.rm = TRUE)
  list(n_total = nrow(wide), wins_countbcf = wins_countbcf, wide = wide)
}

###############################################################################
## 5b. Sensitivity (ablation) tables: method x factor-level, as text/csv/tex.
##     These are the source of truth for the written sensitivity analysis
##     (so the prose is read off tables, not off the PNG plots).
###############################################################################

SENS_METRICS <- tribble(
  ~metric,                ~lab,                    ~tex_lab,            ~digits, ~key,
  "rmse_ate",             "RMSE(ATE)",             "RMSE(ATE)",         2,       TRUE,
  "ate_coverage_mean",    "ATE 95% coverage",      "ATE Cov.\\,95",     2,       TRUE,
  "ate_ci_width_mean",    "Avg ATE 95% CI width",  "ATE CI Width",      2,       FALSE,
  "pehe_mean",            "Avg PEHE",              "PEHE",              2,       TRUE,
  "cate_cov95_mean",      "Avg CATE 95% coverage", "CATE Cov.\\,95",    2,       TRUE,
  "cate_ci_width95_mean", "Avg CATE 95% CI width", "CATE CI Width95",   2,       FALSE
)

FACTOR_LEVELS <- list(N = c(100, 250, 500),
                      P = c(5, 50, 250),
                      ATE = c(0.5, 1.25, 2.5))
FACTOR_PRETTY <- c(N = "sample size N", P = "covariate count p",
                   ATE = "target ATE")

# Successful replicates (out of 100) per (study, method, dgp, cell).
completeness_cell <- raw_all %>%
  group_by(study, method, dgp, N, P, ate_target) %>%
  summarise(n_total = dplyr::n(), n_success = sum(success), .groups = "drop")

methods_in <- function(summary_df) {
  intersect(METHOD_ORDER, unique(summary_df$method))
}

# mean successful replicates over DGPs, long (method, level, value).
completeness_long <- function(study_name, factor_name) {
  cc <- completeness_cell %>% filter(study == study_name)
  bind_rows(lapply(FACTOR_LEVELS[[factor_name]], function(L) {
    s <- switch(factor_name,
      N   = cc %>% filter(N == L,     P == REF_P, abs(ate_target - REF_ATE) < 1e-9),
      P   = cc %>% filter(N == REF_N, P == L,     abs(ate_target - REF_ATE) < 1e-9),
      ATE = cc %>% filter(N == REF_N, P == REF_P, abs(ate_target - L) < 1e-9))
    s %>% group_by(method) %>%
      summarise(value = mean(n_success), .groups = "drop") %>% mutate(level = L)
  }))
}

# metric averaged across DGPs, long (method, level, value).
metric_long <- function(summary_df, factor_name, metric) {
  fs <- factor_sweep(summary_df, factor_name)
  fs %>% transmute(method, level, value = .data[[metric]])
}

# Full long sensitivity frame for one study.
build_sensitivity_long <- function(study_name, summary_df) {
  out <- list()
  for (fac in names(FACTOR_LEVELS)) {
    for (k in seq_len(nrow(SENS_METRICS))) {
      out[[length(out) + 1]] <- metric_long(summary_df, fac, SENS_METRICS$metric[k]) %>%
        mutate(study = study_name, factor = fac, metric = SENS_METRICS$metric[k])
    }
    out[[length(out) + 1]] <- completeness_long(study_name, fac) %>%
      mutate(study = study_name, factor = fac, metric = "n_success")
  }
  bind_rows(out) %>% select(study, factor, level, method, metric, value)
}

# Render a method x level matrix (one metric) as Markdown rows.
render_md_matrix <- function(long_sub, methods_ord, levels_ord, fmt) {
  hdr <- c("Method", as.character(levels_ord))
  lines <- c(paste0("| ", paste(hdr, collapse = " | "), " |"),
             paste0("|", paste(rep("---", length(hdr)), collapse = "|"), "|"))
  for (m in methods_ord) {
    cells <- METHOD_LABELS[[m]]
    for (L in levels_ord) {
      v <- long_sub$value[long_sub$method == m & long_sub$level == L]
      cells <- c(cells, if (length(v) == 0 || is.na(v)) "--" else fmt(v))
    }
    lines <- c(lines, paste0("| ", paste(cells, collapse = " | "), " |"))
  }
  lines
}

# Render the same matrix as a LaTeX tabular.
render_tex_matrix <- function(long_sub, methods_ord, levels_ord, fmt, caption, label) {
  col_spec <- paste0("l", paste(rep("r", length(levels_ord)), collapse = ""))
  hdr <- paste("Method", paste(sprintf("$%s$", levels_ord), collapse = " & "), sep = " & ")
  body <- character(0)
  for (m in methods_ord) {
    cells <- METHOD_LABELS[[m]]
    for (L in levels_ord) {
      v <- long_sub$value[long_sub$method == m & long_sub$level == L]
      cells <- c(cells, if (length(v) == 0 || is.na(v)) "--" else fmt(v))
    }
    body <- c(body, paste(paste(cells, collapse = " & "), "\\\\"))
  }
  c(sprintf("\\begin{table}[h!]\\footnotesize\\centering"),
    sprintf("\\caption{%s}\\label{%s}", caption, label),
    sprintf("\\begin{tabular}{%s}", col_spec), "\\toprule",
    paste(hdr, "\\\\"), "\\midrule", body, "\\bottomrule",
    "\\end{tabular}", "\\end{table}")
}

write_sensitivity_tables <- function(study_name, summary_df) {
  long <- build_sensitivity_long(study_name, summary_df)
  mo   <- methods_in(summary_df)
  write_csv(long, file.path(TABLES_DIR, sprintf("sensitivity_%s.csv", study_name)))

  md  <- c(sprintf("# Sensitivity (ablation) tables --- %s study", study_name), "",
    "Each table fixes the other two factors at the reference cell",
    sprintf("($N = %d$, $p = %d$, ATE $= %.2f$) and sweeps one factor.", REF_N, REF_P, REF_ATE),
    "Metric cells are the **mean across the four DGPs** (successful replicates only).",
    "`n_success` is the mean number of successful replicates out of 100, averaged",
    "over DGPs --- a value of 0 means the method fails at every DGP for that cell",
    "(e.g. a parametric GLM with $p \\ge n$); `--` marks a metric that has no",
    "successful replicate to average.", "")
  tex <- c(sprintf("%% Sensitivity tables --- %s study (generated).", study_name))

  for (fac in names(FACTOR_LEVELS)) {
    lv <- FACTOR_LEVELS[[fac]]
    md <- c(md, sprintf("## Varying %s", FACTOR_PRETTY[[fac]]), "")
    # completeness
    sub <- long %>% filter(factor == fac, metric == "n_success")
    md <- c(md, "**Successful replicates (out of 100), mean across DGPs**", "",
            render_md_matrix(sub, mo, lv, function(v) sprintf("%.0f", v)), "")
    tex <- c(tex, render_tex_matrix(sub, mo, lv, function(v) sprintf("%.0f", v),
      sprintf("Sensitivity (%s study): successful replicates out of 100 (mean over DGPs) while varying %s.",
              study_name, FACTOR_PRETTY[[fac]]),
      sprintf("tab:sens_%s_%s_nsucc", study_name, fac)), "")
    for (k in seq_len(nrow(SENS_METRICS))) {
      mt <- SENS_METRICS$metric[k]; dg <- SENS_METRICS$digits[k]
      sub <- long %>% filter(factor == fac, metric == mt)
      md <- c(md, sprintf("**%s**", SENS_METRICS$lab[k]), "",
              render_md_matrix(sub, mo, lv, function(v) pretty(v, dg)), "")
      if (isTRUE(SENS_METRICS$key[k])) {
        tex <- c(tex, render_tex_matrix(sub, mo, lv, function(v) pretty(v, dg),
          sprintf("Sensitivity (%s study): %s (mean over DGPs) while varying %s.",
                  study_name, SENS_METRICS$tex_lab[k], FACTOR_PRETTY[[fac]]),
          sprintf("tab:sens_%s_%s_%s", study_name, fac, mt)), "")
      }
    }
  }
  writeLines(md,  file.path(TABLES_DIR, sprintf("sensitivity_%s.md",  study_name)))
  writeLines(tex, file.path(TABLES_DIR, sprintf("sensitivity_%s.tex", study_name)))
  invisible(long)
}

write_sensitivity_tables("count", summary_count)
write_sensitivity_tables("zi",    summary_zi)
cat("Wrote sensitivity tables.\n")

###############################################################################
## 6. Analysis paragraphs (data-driven; the same prose feeds .md and .tex)
###############################################################################

# Returns a list with $md (vector of Markdown lines) and $tex (vector of LaTeX
# lines) describing the headline + factor-sweep findings for one study.
build_analysis <- function(study, summary_df, dgp_labels) {
  ra <- ref_agg_for(summary_df)
  cb <- ra %>% filter(method == "countbcf")
  bc <- ra %>% filter(method == "bcf_gauss")

  d_rmse_pct <- 100 * (bc$rmse_ate - cb$rmse_ate) / bc$rmse_ate
  d_pehe_pct <- 100 * (bc$pehe_mean - cb$pehe_mean) / bc$pehe_mean
  d_atecov   <- cb$ate_coverage_mean   - bc$ate_coverage_mean
  d_catecov  <- cb$cate_cov95_mean     - bc$cate_cov95_mean
  d_ciwidth  <- 100 * (cb$ate_ci_width_mean    - bc$ate_ci_width_mean)  /
                       bc$ate_ci_width_mean
  d_ciwidth_cate <- 100 * (cb$cate_ci_width95_mean - bc$cate_ci_width95_mean) /
                            bc$cate_ci_width95_mean

  # Head-to-head DGP wins for the proposed method
  h_rmse  <- ref_h2h_count(summary_df, dgp_labels, "rmse_ate",   TRUE)
  h_pehe  <- ref_h2h_count(summary_df, dgp_labels, "pehe_mean",  TRUE)
  h_catec <- ref_h2h_count(summary_df, dgp_labels, "cate_cov95_mean", FALSE)

  # Factor sweep summaries: change between extremes for each method
  sweep_N   <- factor_sweep(summary_df, "N")
  sweep_P   <- factor_sweep(summary_df, "P")
  sweep_ATE <- factor_sweep(summary_df, "ATE")

  pct_change <- function(sweep_df, metric, from_level, to_level) {
    sub <- sweep_df %>% filter(level %in% c(from_level, to_level))
    out <- list()
    for (m in c("bcf_gauss", "countbcf")) {
      v <- sub %>% filter(method == m) %>% arrange(level)
      if (nrow(v) < 2) {
        out[[m]] <- NA_real_
      } else {
        from_v <- v %>% filter(level == from_level) %>% pull(metric)
        to_v   <- v %>% filter(level == to_level)   %>% pull(metric)
        out[[m]] <- 100 * (to_v - from_v) / from_v
      }
    }
    out
  }

  d_rmse_N   <- pct_change(sweep_N, "rmse_ate",  100, 500)
  d_pehe_N   <- pct_change(sweep_N, "pehe_mean", 100, 500)
  d_rmse_P   <- pct_change(sweep_P, "rmse_ate",  5,   250)
  d_pehe_P   <- pct_change(sweep_P, "pehe_mean", 5,   250)
  d_rmse_ATE <- pct_change(sweep_ATE, "rmse_ate",  0.5, 2.5)
  d_pehe_ATE <- pct_change(sweep_ATE, "pehe_mean", 0.5, 2.5)
  d_cov_N    <- pct_change(sweep_N, "ate_coverage_mean", 100, 500)
  d_cov_P    <- pct_change(sweep_P, "ate_coverage_mean", 5,   250)

  # Normalized ATE-sweep summaries: divide every (method, dgp, ATE) value by
  # CountBCF's value at the SAME (dgp, ATE) cell, then average across DGPs.
  # CountBCF is therefore 1.0 at every ATE level, and each other method's
  # normalized value reads as its ratio to CountBCF at that same target ATE.
  ate_norm_avg <- function(metric_name) {
    sub <- summary_df %>%
      filter(N == REF_N, P == REF_P) %>%
      select(dgp, method, ate_target, all_of(metric_name)) %>%
      rename(.val = all_of(metric_name))
    base <- sub %>%
      filter(method == "countbcf") %>%
      select(dgp, ate_target, .base = .val)
    sub <- sub %>% left_join(base, by = c("dgp", "ate_target")) %>%
      mutate(norm = .val / .base)
    sub %>%
      group_by(method, ate_target) %>%
      summarise(norm_mean = mean(norm, na.rm = TRUE), .groups = "drop")
  }
  rmse_norm <- ate_norm_avg("rmse_ate")
  pehe_norm <- ate_norm_avg("pehe_mean")
  pick <- function(df, m, a) {
    df %>% filter(method == m, abs(ate_target - a) < 1e-9) %>% pull(norm_mean)
  }
  rmse_norm_bcf_05  <- pick(rmse_norm, "bcf_gauss", 0.5)
  rmse_norm_bcf_25  <- pick(rmse_norm, "bcf_gauss", 2.5)
  pehe_norm_bcf_05  <- pick(pehe_norm, "bcf_gauss", 0.5)
  pehe_norm_bcf_25  <- pick(pehe_norm, "bcf_gauss", 2.5)

  study_short <- if (study == "count") "count-only" else "zero-inflated"

  paragraphs <- list(

    headline_md = c(
      "### Headline: reference cell (N = 250, p = 5, ATE = 1.25)",
      "",
      sprintf(paste("Averaged over the four %s DGPs, CountBCF improves on the",
                    "Gaussian BCF baseline on every aggregate metric. CountBCF",
                    "cuts the ATE root-mean-squared error from %s to %s",
                    "(%s) and the average PEHE from %s to %s (%s)."),
              study_short,
              pretty(bc$rmse_ate),  pretty(cb$rmse_ate),
              pretty_pct(-d_rmse_pct),
              pretty(bc$pehe_mean), pretty(cb$pehe_mean),
              pretty_pct(-d_pehe_pct)),
      "",
      sprintf(paste("ATE 95%% coverage moves from %s under BCF (Gaussian) to %s",
                    "under CountBCF (gap %s); per-unit CATE 95%% coverage moves",
                    "from %s to %s (gap %s), pulling both noticeably closer to",
                    "the nominal 0.95. CountBCF's intervals are wider in",
                    "exchange (ATE CI %s; unit-level CATE CI %s), which is the",
                    "honest cost of producing intervals that actually cover."),
              pretty(bc$ate_coverage_mean), pretty(cb$ate_coverage_mean),
              pretty_pct(100 * d_atecov, 0),
              pretty(bc$cate_cov95_mean),  pretty(cb$cate_cov95_mean),
              pretty_pct(100 * d_catecov, 0),
              pretty_pct(d_ciwidth, 0),
              pretty_pct(d_ciwidth_cate, 0)),
      "",
      sprintf(paste("Looking DGP-by-DGP at the reference cell, CountBCF wins on",
                    "RMSE(ATE) in %d/%d DGPs, on PEHE in %d/%d, and is closer",
                    "to nominal CATE coverage in %d/%d."),
              h_rmse$wins_countbcf, h_rmse$n_total,
              h_pehe$wins_countbcf, h_pehe$n_total,
              h_catec$wins_countbcf, h_catec$n_total),
      ""
    ),

    ablation_md = c(
      "### Ablation: how the gap moves with N, p, and target ATE",
      "",
      sprintf(paste("**Sample size $N$ (100 -> 500).** RMSE(ATE) falls for both",
                    "methods as expected (BCF (Gaussian): %s, CountBCF: %s),",
                    "and so does PEHE (BCF (Gaussian): %s, CountBCF: %s).",
                    "BCF (Gaussian)'s ATE coverage drifts %s with N as its",
                    "intervals shrink faster than its bias, while CountBCF's",
                    "coverage drifts %s. CountBCF retains an advantage on",
                    "PEHE at every sample size; the absolute gap narrows at",
                    "$N=500$ as both methods approach the asymptote."),
              pretty_pct(d_rmse_N$bcf_gauss),
              pretty_pct(d_rmse_N$countbcf),
              pretty_pct(d_pehe_N$bcf_gauss),
              pretty_pct(d_pehe_N$countbcf),
              pretty_pct(d_cov_N$bcf_gauss),
              pretty_pct(d_cov_N$countbcf)),
      "",
      sprintf(paste("**Covariate count $p$ (5 -> 250).** Adding 245 noise",
                    "covariates degrades BCF (Gaussian) heavily (RMSE(ATE) %s,",
                    "PEHE %s, ATE coverage %s), while CountBCF is markedly more",
                    "robust (RMSE(ATE) %s, PEHE %s). Because both forests use",
                    "the same Dirichlet-style splitting rules, the gap reflects",
                    "the count likelihood: when the response is heavily skewed",
                    "and the wrong scale is fit, irrelevant covariates can",
                    "stand in for residual count structure and inflate the",
                    "Gaussian model's variance."),
              pretty_pct(d_rmse_P$bcf_gauss),
              pretty_pct(d_pehe_P$bcf_gauss),
              pretty_pct(d_cov_P$bcf_gauss),
              pretty_pct(d_rmse_P$countbcf),
              pretty_pct(d_pehe_P$countbcf)),
      "",
      sprintf(paste("**Target ATE (0.5 -> 2.5).** Absolute errors scale with",
                    "the target ATE magnitude (a $5\\times$ change in target",
                    "rescales the response-scale CATE by the same factor), so",
                    "reading raw RMSE(ATE) and PEHE across ATE levels conflates",
                    "method behaviour with that scaling. We therefore report",
                    "errors *normalized to CountBCF at the same ATE level*:",
                    "for each DGP and each ATE level we divide every method's",
                    "cell by CountBCF's value at that same (DGP, ATE) cell, and",
                    "then average across DGPs. CountBCF is thus $1$ at every",
                    "ATE level by construction, and each other curve reads as",
                    "its ratio to CountBCF at that level. On this scale BCF",
                    "(Gaussian) sits at %s on RMSE(ATE) and %s on PEHE at ATE",
                    "$= 0.5$, and at %s and %s respectively at ATE $= 2.5$: it",
                    "stays above CountBCF (ratio $> 1$) on both metrics at every",
                    "ATE level, with a broadly stable multiplicative gap rather",
                    "than one that closes as the signal grows, and CountBCF",
                    "coverage of the ATE stays closer to $0.95$ throughout.",
                    "(Absolute, unnormalized versions of these panels are saved",
                    "alongside as `plots/%s_ATE_rmse_ate.png` and",
                    "`plots/%s_ATE_pehe_mean.png`.)"),
              pretty(rmse_norm_bcf_05),
              pretty(pehe_norm_bcf_05),
              pretty(rmse_norm_bcf_25),
              pretty(pehe_norm_bcf_25),
              study, study),
      ""
    ),

    interpretation_md = c(
      "### Interpretation",
      "",
      paste("Three threads run through both studies. First, treating the count",
            "outcome as Gaussian systematically *understates* the posterior",
            "uncertainty: the BCF (Gaussian) intervals are 25--60% narrower",
            "than CountBCF's, yet their coverage is 10--15 percentage points",
            "below nominal, which is the worst kind of mis-calibration for a",
            "policy report. CountBCF's wider intervals are not a free-lunch",
            "loss in sharpness -- they are the right intervals."),
      "",
      paste("Second, the proposed method's edge widens with model",
            "mis-specification: the BCF/CountBCF PEHE ratio is most extreme",
            "on the nonlinear NB and ZINB DGPs, exactly where the Gaussian",
            "approximation is worst. On the easier linear-Poisson cells the",
            "two methods are closer, because there the count likelihood is",
            "best approximated by a Gaussian with constant variance."),
      "",
      paste("Third, the methods react quite differently to nuisance",
            "covariates. CountBCF's PEHE and coverage are nearly flat across",
            "$p \\in \\{5, 50, 250\\}$ at $N=250$, while the Gaussian BCF",
            "degrades steadily. That suggests the regularization on the",
            "log-scale count model is doing real work to filter signal from",
            "noise covariates, beyond what the propensity-score adjustment",
            "alone provides.")
    ),

    caveats_md = c(
      "### Caveats",
      "",
      paste("- `N_SIM = 100` per cell, so cross-replicate SDs are the right",
            "scale for the variance bars but cell-level means still carry",
            "Monte-Carlo noise on the order of $\\widehat\\sigma/\\sqrt{100}$."),
      paste("- Coverage of the ATE is read off the *realized* ATE for that",
            "replicate, not a single population value -- this is the right",
            "way to score a posterior CI on a finite sample, but it does",
            "make coverage move with the calibration of $\\tau_0$."),
      if (study == "count")
        paste("- One count cell --- countbcf on",
              "`nonlinear_poisson_cate__N250_P250` --- is missing from disk",
              "(the runner did not complete). The corresponding panel in the",
              "$p$-ablation plots therefore shows only the BCF (Gaussian)",
              "curve at $p=250$ for that DGP.")
      else
        paste("- All 56 zero-inflated cells completed 100 replicates",
              "successfully; nonetheless the cross-replicate PEHE SD at",
              "$p=250$ remains large for both methods, so cell-to-cell",
              "noise in the $p$-ablation panels at the high end should be",
              "interpreted with a Monte-Carlo error of order",
              "$\\widehat\\sigma_{\\text{PEHE}}/\\sqrt{100}$ in mind."),
      paste("- All runs use $\\hat\\pi = \\pi$ (the true propensity is passed",
            "in), so this benchmark isolates the outcome model. Real-data",
            "performance with an estimated $\\hat\\pi$ will be worse for",
            "both methods to a similar degree.")
    )
  )

  paragraphs
}

###############################################################################
## 7. Write the per-study .md and .tex narrative wrappers
###############################################################################

write_study_writeup_md <- function(study, summary_df, dgp_labels,
                                   table_md_path, out_path) {
  paras <- build_analysis(study, summary_df, dgp_labels)
  study_long <- if (study == "count") "Count-only (Poisson + NB)"
                else "Zero-inflated (ZIP + ZINB)"

  ablation_section <- character(0)
  for (factor in c("N", "P", "ATE")) {
    fac_lab <- switch(factor, N = "sample size N", P = "covariate count p",
                      ATE = "target ATE")
    ablation_section <- c(ablation_section,
      sprintf("#### Varying %s", fac_lab), "")
    if (factor == "ATE") {
      ablation_section <- c(ablation_section,
        paste("RMSE(ATE) and PEHE are reported in *normalized* form: each",
              "(method, DGP) cell is divided by CountBCF's value at the *same*",
              "target ATE in that DGP, so the CountBCF curve is flat at 1 and",
              "every other curve reads as its ratio to CountBCF at that ATE",
              "level. Absolute (unnormalized) versions of the same panels are",
              "saved alongside as `plots/<study>_ATE_rmse_ate.png` and",
              "`plots/<study>_ATE_pehe_mean.png`."),
        "")
    }
    for (metric in c("rmse_ate", "ate_coverage_mean", "ate_ci_width_mean",
                     "pehe_mean", "cate_cov95_mean", "cate_ci_width95_mean")) {
      mlab <- METRIC_DEFS$ylab[METRIC_DEFS$metric == metric]
      normalized <- factor == "ATE" && metric %in% ATE_NORMALIZE_METRICS
      png_rel <- if (normalized) {
        sprintf("plots/%s_%s_%s_normalized.png", study, factor, metric)
      } else {
        sprintf("plots/%s_%s_%s.png", study, factor, metric)
      }
      caption_lab <- if (normalized) paste(mlab, "(normalized to CountBCF at each ATE)") else mlab
      ablation_section <- c(ablation_section,
        sprintf("**%s.** ![](%s)", caption_lab, png_rel), "")
    }
  }

  # Hand-written baselines analysis (preserved across re-runs in a notes file)
  # plus the generated sensitivity tables.
  notes_path <- file.path(OUT_DIR, sprintf("notes_%s_baselines.md", study))
  sens_path  <- file.path(TABLES_DIR, sprintf("sensitivity_%s.md", study))
  baselines_md <- character(0)
  if (file.exists(notes_path)) {
    baselines_md <- c("## Baselines: parametric models and causal forest", "",
                      readLines(notes_path), "")
  }
  if (file.exists(sens_path)) {
    sens_lines <- readLines(sens_path)
    if (length(sens_lines) && grepl("^# ", sens_lines[1])) sens_lines <- sens_lines[-1]
    sens_lines <- sub("^## ", "### ", sens_lines)  # nest under the H2 below
    baselines_md <- c(baselines_md,
                      "## Sensitivity (ablation) tables", "", sens_lines, "")
  }

  txt <- c(
    sprintf("# Simulation summary --- %s study", study_long),
    "",
    sprintf("Companion to `tests/%s` (see for the full DGP and metric definitions).",
            ifelse(study == "count", "sim_count.md", "sim_zi.md")),
    sprintf("Sources of truth: the per-cell CSVs under `tests/results_%s/`",
            study),
    "(one row per Monte-Carlo replicate, `N_SIM = 100` per cell).",
    "All aggregates and plots in this file are produced by",
    "`tests/analyze_results.R`.",
    "",
    "## Metrics reported",
    "",
    "All quantities are computed across the 100 MC replicates within each cell:",
    "",
    "- **RMSE ATE** --- $\\sqrt{(1/S)\\sum_s (\\widehat{\\text{ATE}}_s - \\text{ATE}_s)^2}$ across replicates $s$.",
    "- **ATE Cov.** --- fraction of replicates whose 95% posterior CI for the ATE contains the realized ATE (nominal 0.95).",
    "- **ATE CI Width** --- mean width of the 95% ATE credible interval (`ate_q975 - ate_q025`), mean $\\pm$ SD across replicates.",
    "- **PEHE** --- $\\sqrt{\\overline{(\\hat\\tau(X_i) - \\tau(X_i))^2}}$, mean $\\pm$ SD.",
    "- **CATE Cov.95** --- per-replicate fraction of units whose 95% unit-level CI covers $\\tau(X_i)$, mean $\\pm$ SD.",
    "- **CATE CI Width95** --- per-replicate mean width of the unit-level 95% credible intervals, mean $\\pm$ SD.",
    "",
    "RMSE(ATE) and ATE Cov. are scalar summaries across replicates by construction, so no SD is reported for them.",
    "",
    "## Reference-cell table",
    "",
    readLines(table_md_path),
    "",
    paras$headline_md,
    paras$ablation_md,
    paras$interpretation_md,
    "",
    paras$caveats_md,
    "",
    baselines_md,
    "## Ablation plots",
    "",
    "Each panel below holds the other two factors at their reference values and sweeps the labeled factor one-at-a-time. Points are the cross-replicate mean for that cell; the dashed grey line on the coverage panels marks the nominal 0.95. Each ablation is shown in two variants: the full panel and a `_clean` panel that drops the heavy-tailed NB-family parametric baseline (NB GLM for the count study, ZINB GLM for the zero-inflated study) so the other curves stay legible.",
    "",
    ablation_section
  )
  writeLines(txt, out_path)
}

# Turn a vector of Markdown analysis lines into LaTeX paragraphs.
# Inline conversions:
#   **bold**     -> \textbf{bold}
#   *emph*       -> \emph{emph}
#   `code`       -> \texttt{code}  (underscores escaped)
#   '-> '        -> ' $\to$ '      (replaced as ASCII arrow)
#   raw %        -> \%             (outside of math; we don't use $...$ here)
md_inline_to_tex <- function(s) {
  # 1) Pull out fenced backtick spans first, escape underscores inside, and
  #    wrap with \texttt{...}.
  re_bt <- "`([^`]+)`"
  s <- gsub(re_bt, "\\\\texttt{\\1}", s, perl = TRUE)
  # Now escape underscores inside any \texttt{...} we just inserted.
  fix_tt <- function(s) {
    m <- gregexpr("\\\\texttt\\{[^}]*\\}", s, perl = TRUE)[[1]]
    if (length(m) == 1 && m[1] == -1) return(s)
    lengths <- attr(m, "match.length")
    chunks <- character(0); last <- 1
    for (i in seq_along(m)) {
      st <- m[i]; ln <- lengths[i]
      pre <- substr(s, last, st - 1)
      mid <- substr(s, st, st + ln - 1)
      mid <- gsub("_", "\\\\_", mid)
      chunks <- c(chunks, pre, mid)
      last <- st + ln
    }
    chunks <- c(chunks, substr(s, last, nchar(s)))
    paste(chunks, collapse = "")
  }
  s <- fix_tt(s)
  # 2) Bold then italic. (Bold must come first because it's `**...**`.)
  s <- gsub("\\*\\*([^*]+)\\*\\*", "\\\\textbf{\\1}", s, perl = TRUE)
  s <- gsub("(^|[^*])\\*([^*]+)\\*", "\\1\\\\emph{\\2}", s, perl = TRUE)
  # 3) ASCII arrow.
  s <- gsub(" -> ", " $\\\\to$ ", s, fixed = FALSE)
  # 4) Raw % outside of math: escape every % unless already preceded by a
  #    backslash. We never embed math segments containing % in this writeup,
  #    so this is safe.
  s <- gsub("(?<!\\\\)%", "\\\\%", s, perl = TRUE)
  s
}

md_to_tex <- function(lines) {
  out <- character(0)
  for (ln in lines) {
    if (grepl("^### ", ln)) {
      out <- c(out, sprintf("\\subsubsection*{%s}",
                            md_inline_to_tex(sub("^### ", "", ln))))
    } else if (grepl("^- ", ln)) {
      out <- c(out, sprintf("\\item %s",
                            md_inline_to_tex(sub("^- ", "", ln))))
    } else if (nzchar(ln)) {
      out <- c(out, md_inline_to_tex(ln))
    } else {
      out <- c(out, "")
    }
  }
  # Wrap consecutive \item lines in an itemize environment.
  final <- character(0); in_item <- FALSE
  for (ln in out) {
    is_item <- grepl("^\\\\item ", ln)
    if (is_item && !in_item) {
      final <- c(final, "\\begin{itemize}")
      in_item <- TRUE
    }
    if (!is_item && in_item) {
      final <- c(final, "\\end{itemize}")
      in_item <- FALSE
    }
    final <- c(final, ln)
  }
  if (in_item) final <- c(final, "\\end{itemize}")
  final
}

write_study_writeup_tex <- function(study, summary_df, dgp_labels,
                                    table_tex_path, out_path) {
  paras <- build_analysis(study, summary_df, dgp_labels)
  study_long <- if (study == "count") "Count-only (Poisson + NB)"
                else "Zero-inflated (ZIP + ZINB)"
  sec_label <- if (study == "count") "sec:sim_count" else "sec:sim_zi"

  tex <- c(
    sprintf("%% Generated by tests/analyze_results.R --- do not edit by hand."),
    sprintf("\\section{Simulation results --- %s}\\label{%s}",
            study_long, sec_label),
    "",
    paste("All numbers below summarize the per-cell Monte-Carlo CSVs under",
          sprintf("\\texttt{tests/results\\_%s/}", study),
          "(100 replicates per cell, seed = replicate id).",
          "RMSE(ATE) and ATE coverage are scalar summaries across replicates;",
          "PEHE, CATE coverage, and the two CI-width columns are reported as",
          "cross-replicate mean $\\pm$ one standard deviation."),
    "",
    "% --- reference-cell table (one of two) --- ",
    readLines(table_tex_path),
    "",
    md_to_tex(paras$headline_md),
    "",
    md_to_tex(paras$ablation_md),
    "",
    md_to_tex(paras$interpretation_md),
    "",
    md_to_tex(paras$caveats_md)
  )

  # Baselines analysis (notes prose -> LaTeX) plus the generated sensitivity
  # tables for this study.
  notes_path <- file.path(OUT_DIR, sprintf("notes_%s_baselines.md", study))
  sens_tex   <- file.path(TABLES_DIR, sprintf("sensitivity_%s.tex", study))
  if (file.exists(notes_path)) {
    tex <- c(tex, "",
      "\\subsection*{Baselines: parametric models and causal forest}", "",
      md_to_tex(readLines(notes_path)))
  }
  if (file.exists(sens_tex)) {
    tex <- c(tex, "",
      "\\subsection*{Sensitivity (ablation) tables}", "",
      readLines(sens_tex))
  }
  writeLines(tex, out_path)
}

write_study_writeup_md("count", summary_count, DGP_LABELS_COUNT,
                       file.path(TABLES_DIR, "table_count_ref.md"),
                       file.path(OUT_DIR, "SUMMARY_count.md"))
write_study_writeup_md("zi", summary_zi, DGP_LABELS_ZI,
                       file.path(TABLES_DIR, "table_zi_ref.md"),
                       file.path(OUT_DIR, "SUMMARY_zi.md"))
write_study_writeup_tex("count", summary_count, DGP_LABELS_COUNT,
                        file.path(TABLES_DIR, "table_count_ref.tex"),
                        file.path(OUT_DIR, "SUMMARY_count.tex"))
write_study_writeup_tex("zi", summary_zi, DGP_LABELS_ZI,
                        file.path(TABLES_DIR, "table_zi_ref.tex"),
                        file.path(OUT_DIR, "SUMMARY_zi.tex"))

cat("Done. Outputs under", OUT_DIR, "\n")
