# =============================================================================
# Question 3a: TLG - Adverse Events Summary Table using {gtsummary}
# =============================================================================
# Objective: Create a regulatory-compliant summary table of Treatment-Emergent
#            Adverse Events (TEAEs) using {gtsummary}, similar to FDA Table 10.
#
# Input:     pharmaverseadam::adae, pharmaverseadam::adsl
# Output:    ae_summary_table.html
#
# Table structure:
#   - Rows:    Primary System Organ Class (AESOC) with nested AETERM
#   - Columns: Treatment groups (ACTARM) + Total
#   - Values:  Count (n) and percentage (%) of subjects
#   - Sort:    Descending frequency within each SOC
#
# Author:    Mridul K. Thomas
# Date:      03-2026
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load Required Libraries
# -----------------------------------------------------------------------------
library(gtsummary)       # Clinical table creation
library(pharmaverseadam) # ADaM example datasets
library(dplyr)           # Data manipulation
library(tidyr)           # Data reshaping
library(gt)              # Table rendering (used by gtsummary)
library(stringr)         # String operations
library(tidytable)
library(tibble)
library(purrr)

# -----------------------------------------------------------------------------
# 2. Load and Inspect Data
# -----------------------------------------------------------------------------
adae <- pharmaverseadam::adae
adsl <- pharmaverseadam::adsl

cat("=== ADAE Dimensions:", nrow(adae), "rows x", ncol(adae), "cols ===\n")
cat("=== ADSL Dimensions:", nrow(adsl), "rows x", ncol(adsl), "cols ===\n")

# Inspect treatment emergent flag
cat("\n=== TRTEMFL Distribution ===\n")
print(table(adae$TRTEMFL, useNA = "always"))

cat("\n=== ACTARM Distribution in ADSL ===\n")
print(table(adsl$ACTARM, useNA = "always"))

# -----------------------------------------------------------------------------
# 3. Prepare Data for Table
# -----------------------------------------------------------------------------
# Filter to Treatment-Emergent AEs only
adae_te <- adae %>%
  filter(TRTEMFL == "Y")

cat("\n=== TEAE count:", nrow(adae_te), "records ===\n")

# Get denominator: number of subjects per treatment arm from ADSL
# (use all subjects, not just those with AEs, for correct % calculation)
n_by_arm <- adsl %>%
  count(ACTARM, name = "N_ARM")

n_total <- nrow(adsl)

cat("\n=== Subjects per arm ===\n")
print(n_by_arm)

# Build subject-level AE indicator dataset:
# One row per subject x SOC x term, flagged Y/N for AE occurrence
# This ensures percentages are of subjects (not events)
subj_ae <- adae_te %>%
  # Get distinct subject x SOC x term combinations (count subjects, not events)
  distinct(USUBJID, AESOC, AETERM) %>%
  mutate(has_ae = 1L)

# Build the full subject x term grid to include zeros
all_subjects <- adsl %>%
  select(USUBJID, ACTARM)

all_terms <- subj_ae %>%
  distinct(AESOC, AETERM)

# Create full grid: every subject x every (SOC, term) combination
full_grid <- crossing(all_subjects, all_terms) %>%
  left_join(subj_ae, by = c("USUBJID", "AESOC", "AETERM")) %>%
  mutate(has_ae = replace_na(has_ae, 0L))

# -----------------------------------------------------------------------------
# 4. Compute Summary Statistics
# -----------------------------------------------------------------------------
# Count subjects with each AE by arm and overall
ae_summary <- full_grid %>%
  group_by(AESOC, AETERM, ACTARM) %>%
  summarise(
    n_subj = sum(has_ae),
    .groups = "drop"
  ) %>%
  # Add arm denominators
  left_join(n_by_arm, by = "ACTARM") %>%
  mutate(pct = round(100 * n_subj / N_ARM, 1))

# Total column
ae_total <- full_grid %>%
  group_by(AESOC, AETERM) %>%
  summarise(
    n_subj_total = sum(has_ae),
    .groups = "drop"
  ) %>%
  mutate(
    N_ARM_total = n_total,
    pct_total   = round(100 * n_subj_total / N_ARM_total, 1)
  )

# Combine
ae_combined <- ae_summary %>%
  pivot_wider(
    id_cols    = c(AESOC, AETERM),
    names_from = ACTARM,
    values_from = c(n_subj, pct)
  ) %>%
  left_join(ae_total, by = c("AESOC", "AETERM"))

# Sort: by SOC frequency descending, then by AETERM frequency descending
soc_order <- ae_total %>%
  group_by(AESOC) %>%
  summarise(soc_total = sum(n_subj_total), .groups = "drop") %>%
  arrange(desc(soc_total))

ae_combined <- ae_combined %>%
  left_join(soc_order, by = "AESOC") %>%
  arrange(desc(soc_total), AESOC, desc(n_subj_total)) %>%
  select(-soc_total)

# -----------------------------------------------------------------------------
# 5. Create gtsummary Table
# -----------------------------------------------------------------------------
# Approach: Use tbl_summary on subject-level binary flags per arm

# Reshape to wide format for gtsummary
# Create binary flag columns for each AE term
# Use a simpler approach: manual gt table for full control

arm_levels <- sort(unique(adsl$ACTARM))
arm_ns      <- setNames(n_by_arm$N_ARM, n_by_arm$ACTARM)

# Build display data frame
build_display_row <- function(soc, term, data) {
  row_data <- data %>%
    filter(AESOC == soc, AETERM == term)

  cells <- sapply(arm_levels, function(arm) {
    n   <- row_data[[paste0("n_subj_", arm)]]
    pct <- row_data[[paste0("pct_", arm)]]
    if (is.na(n) || n == 0) return("0 (0%)")
    sprintf("%d (%.1f%%)", n, pct)
  })

  total_n   <- row_data$n_subj_total
  total_pct <- row_data$pct_total
  total_cell <- if (is.na(total_n) || total_n == 0) "0 (0%)"
                else sprintf("%d (%.1f%%)", total_n, total_pct)

  c(SOC = soc, TERM = term, cells, Total = total_cell)
}

# Add overall TEAEs row at top
teae_row <- bind_cols(
  tibble(
    SOC  = "Treatment Emergent AEs",
    TERM = ""
  ),
  as_tibble(t(sapply(arm_levels, function(arm) {
    n <- sum(!is.na(adae_te$USUBJID[adae_te$ACTARM == arm]) &
               !duplicated(adae_te[adae_te$ACTARM == arm, "USUBJID"]))
    # Count unique subjects with any TEAE per arm
    n <- adsl %>%
      filter(ACTARM == arm) %>%
      summarise(
        n = sum(USUBJID %in% adae_te$USUBJID)
      ) %>%
      pull(n)
    pct <- round(100 * n / arm_ns[[arm]], 1)
    sprintf("%d\n(%d%%)", n, round(pct))
  })))
) %>%
  mutate(Total = {
    n <- n_distinct(adae_te$USUBJID)
    pct <- round(100 * n / n_total)
    sprintf("%d\n(%d%%)", n, pct)
  })

names(teae_row)[3:(2 + length(arm_levels))] <- arm_levels

# Build all rows
all_rows <- map2_dfr(
  ae_combined$AESOC,
  ae_combined$AETERM,
  ~ as_tibble_row(build_display_row(.x, .y, ae_combined))
)

# -----------------------------------------------------------------------------
# 6. Build GT Table
# -----------------------------------------------------------------------------
# Create header with N per arm
arm_headers <- c(
  "Primary System Organ Class\nReported Term for the Adverse Event",
  sapply(arm_levels, function(arm) {
    sprintf("%s\nN = %d", arm, arm_ns[[arm]])
  }),
  sprintf("Total\nN = %d", n_total)
)

# Combine TEAE header row with detail rows
display_tbl <- bind_rows(teae_row, all_rows)

# Rename columns for gt
names(display_tbl) <- c("SOC", "TERM", arm_levels, "Total")

# Create the GT table
gt_tbl <- display_tbl %>%
  select(-SOC) %>%  # SOC will be shown via row grouping
  gt(rowname_col = "TERM") %>%
  tab_header(
    title    = "Summary of Treatment-Emergent Adverse Events",
    subtitle = "Safety Population"
  ) %>%
  cols_label(
    .list = setNames(
      c(sapply(arm_levels, function(a) sprintf("%s\nN = %d", a, arm_ns[[a]])),
        sprintf("Total\nN = %d", n_total)),
      c(arm_levels, "Total")
    )
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(indent = px(15)),
    locations = cells_stub(rows = TERM != "")
  ) %>%
  tab_footnote(
    footnote = "Percentages are based on number of subjects in the safety population per treatment arm.",
    locations = cells_column_labels(columns = everything())
  ) %>%
  opt_table_font(font = list(google_font("Source Sans Pro"), default_fonts())) %>%
  tab_options(
    table.font.size        = px(12),
    column_labels.font.weight = "bold",
    row.striping.include_table_body = TRUE
  )

# -----------------------------------------------------------------------------
# 7. Alternative: Use gtsummary tbl_summary approach
# -----------------------------------------------------------------------------
# This approach leverages gtsummary's built-in summary mechanisms

# Create long format dataset with binary indicator per subject per AE
ae_for_gtsum <- adsl %>%
  select(USUBJID, ACTARM) %>%
  left_join(
    adae_te %>%
      distinct(USUBJID, AESOC, AETERM) %>%
      mutate(has_ae = TRUE),
    by = "USUBJID"
  ) %>%
  mutate(has_ae = replace_na(has_ae, FALSE))

# Get top AEs by frequency for a cleaner table
top_aes <- ae_for_gtsum %>%
  filter(has_ae) %>%
  count(AETERM, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(AETERM)

ae_for_gtsum_top <- ae_for_gtsum %>%
  filter(is.na(AETERM) | AETERM %in% top_aes) %>%
  mutate(
    AETERM = replace_na(AETERM, "(No AE)"),
    has_ae = as.integer(has_ae)
  )

# Build gtsummary table
tbl_ae <- ae_for_gtsum_top %>%
  filter(AETERM != "(No AE)") %>%
  select(ACTARM, AESOC, AETERM, has_ae) %>%
  tbl_summary(
    by        = ACTARM,
    include   = c(AESOC, AETERM),
    statistic = all_categorical() ~ "{n} ({p}%)",
    label     = list(
      AESOC  ~ "Primary System Organ Class",
      AETERM ~ "Reported Term for the Adverse Event"
    )
  ) %>%
  add_overall(last = TRUE, col_label = "**Total**") %>%
  bold_labels() %>%
  modify_header(
    label ~ "**Primary System Organ Class**\n**Reported Term for the Adverse Event**"
  ) %>%
  modify_caption("**Summary of Treatment-Emergent Adverse Events (Safety Population)**")

# -----------------------------------------------------------------------------
# 8. Save Output
# -----------------------------------------------------------------------------
# Save as HTML
as_gt(tbl_ae) %>%
  gtsave("ae_summary_table.html")

cat("\n✓ AE Summary Table saved to ae_summary_table.html\n")

capture.output(
  {
    cat("Question 3a log\n")
    cat("================\n")
    cat("ADAE rows:", nrow(adae), "\n")
    cat("ADSL rows:", nrow(adsl), "\n")
    cat("TEAE rows:", nrow(adae_te), "\n")
    cat("Total subjects:", n_total, "\n\n")
    
    cat("Subjects per arm:\n")
    print(n_by_arm)
    
    cat("\nTRTEMFL distribution:\n")
    print(table(adae$TRTEMFL, useNA = "always"))
    
    cat("\nACTARM distribution in ADSL:\n")
    print(table(adsl$ACTARM, useNA = "always"))
    
    cat("\nAE summary dimensions:\n")
    cat("Rows:", nrow(ae_combined), " Cols:", ncol(ae_combined), "\n")
    
    cat("\nOutput file:\n")
    cat("ae_summary_table.html\n")
  },
  file = "question_3a_log.txt"
)

cat("\n✓ Log saved to question_3a_log.txt\n")
cat("\n=== Script completed successfully ===\n")
