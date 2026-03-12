# =============================================================================
# Question 2: ADaM ADSL Dataset Creation using {admiral}
# =============================================================================
# Objective: Create an ADSL (Subject-Level Analysis Dataset) using SDTM source
#            data and the {admiral} package from the Pharmaverse ecosystem.
#
# Input datasets:
#   - pharmaversesdtm::dm  (Demographics - basis of ADSL)
#   - pharmaversesdtm::vs  (Vital Signs - for LSTAVLDT)
#   - pharmaversesdtm::ex  (Exposure - for TRTSDTM, TRTEDTM, LSTAVLDT)
#   - pharmaversesdtm::ds  (Disposition - for LSTAVLDT)
#   - pharmaversesdtm::ae  (Adverse Events - for LSTAVLDT)
#
# Custom variables derived:
#   - AGEGR9 / AGEGR9N : Age group (<18, 18-50, >50)
#   - TRTSDTM / TRTSTMF : Treatment start datetime with imputation flag
#   - ITTFL              : Intent-to-Treat flag
#   - LSTAVLDT           : Last known alive date
#
# Author:    Mridul K. Thomas
# Date:      03-2026
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load libraries
# -----------------------------------------------------------------------------
library(admiral)
library(pharmaversesdtm)
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(rlang)

# -----------------------------------------------------------------------------
# 2. Load source datasets
# -----------------------------------------------------------------------------
dm <- pharmaversesdtm::dm %>% convert_blanks_to_na()
vs <- pharmaversesdtm::vs %>% convert_blanks_to_na()
ex <- pharmaversesdtm::ex %>% convert_blanks_to_na()
ds <- pharmaversesdtm::ds %>% convert_blanks_to_na()
ae <- pharmaversesdtm::ae %>% convert_blanks_to_na()

cat("=== Source Data Dimensions ===\n")
cat("DM:", nrow(dm), "rows\n")
cat("VS:", nrow(vs), "rows\n")
cat("EX:", nrow(ex), "rows\n")
cat("DS:", nrow(ds), "rows\n")
cat("AE:", nrow(ae), "rows\n")

# -----------------------------------------------------------------------------
# 3. Initialize ADSL from DM
# -----------------------------------------------------------------------------
adsl <- dm %>%
  select(
    STUDYID, USUBJID, SUBJID, SITEID,
    AGE, AGEU, SEX, RACE, ETHNIC, COUNTRY,
    DMDTC, DMDY,
    ARM, ACTARM, ARMCD, ACTARMCD,
    RFSTDTC, RFENDTC, RFXSTDTC, RFXENDTC, RFICDTC, RFPENDTC,
    DTHDTC, DTHFL,
    any_of("BRTHDTC")
  ) %>%
  mutate(
    RFSTDT = convert_dtc_to_dt(RFSTDTC),
    RFENDT = convert_dtc_to_dt(RFENDTC)
  )

cat("\n=== ADSL initialized from DM ===\n")
cat("Subjects:", nrow(adsl), "\n")

# -----------------------------------------------------------------------------
# 4. Derive AGEGR9 / AGEGR9N
# -----------------------------------------------------------------------------
adsl <- adsl %>%
  mutate(
    AGEGR9 = case_when(
      AGE < 18 ~ "<18",
      AGE >= 18 & AGE <= 50 ~ "18 - 50",
      AGE > 50 ~ ">50",
      TRUE ~ NA_character_
    ),
    AGEGR9N = case_when(
      AGEGR9 == "<18" ~ 1L,
      AGEGR9 == "18 - 50" ~ 2L,
      AGEGR9 == ">50" ~ 3L,
      TRUE ~ NA_integer_
    )
  )

cat("\n=== AGEGR9 Distribution ===\n")
print(table(adsl$AGEGR9, useNA = "always"))

# -----------------------------------------------------------------------------
# 5. Derive ITTFL
# -----------------------------------------------------------------------------
adsl <- adsl %>%
  mutate(
    ITTFL = if_else(!is.na(ARM) & ARM != "", "Y", "N")
  )

cat("\n=== ITTFL Distribution ===\n")
print(table(adsl$ITTFL, useNA = "always"))

# -----------------------------------------------------------------------------
# 6. Prepare EX for treatment date derivations
# -----------------------------------------------------------------------------
# Valid dose:
#   EXDOSE > 0 OR
#   EXDOSE == 0 and EXTRT contains PLACEBO
#
# For TRTSDTM/TRTSTMF the assessment requires a complete datepart of EXSTDTC.
# We therefore precompute EXSTDT and keep only records with complete datepart.
ex_valid <- ex %>%
  filter(
    EXDOSE > 0 |
      (EXDOSE == 0 & str_detect(EXTRT, regex("PLACEBO", ignore_case = TRUE)))
  ) %>%
  mutate(
    EXSTDT = convert_dtc_to_dt(EXSTDTC),
    EXENDT = convert_dtc_to_dt(EXENDTC)
  )

# Start datetime: complete EXSTDTC datepart only, impute missing time only
ex_start <- ex_valid %>%
  filter(!is.na(EXSTDT)) %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST",
    highest_imputation = "h",
    time_imputation = "first",
    flag_imputation = "time"
  )

cat("\n=== Valid EX start records ===\n")
cat("Records:", nrow(ex_start), "\n")

# -----------------------------------------------------------------------------
# 7. Derive TRTSDTM / TRTSTMF
# -----------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_start,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    filter_add = !is.na(EXSTDTM)
  )

cat("\n=== TRTSDTM sample ===\n")
print(head(adsl %>% select(USUBJID, TRTSDTM, TRTSTMF), 5))

# -----------------------------------------------------------------------------
# 8. Derive TRTEDTM
# -----------------------------------------------------------------------------
# Complete EXENDTC datepart only, impute missing time only
ex_end <- ex_valid %>%
  filter(!is.na(EXENDT)) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    highest_imputation = "h",
    time_imputation = "last",
    flag_imputation = "time"
  )

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_end,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(TRTEDTM = EXENDTM),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    filter_add = !is.na(EXENDTM)
  )

# -----------------------------------------------------------------------------
# 9 Derive LSTAVLDT source dates
# -----------------------------------------------------------------------------
# (1) Last complete VS date with valid result
vs_dates <- vs %>%
  mutate(VSDT = convert_dtc_to_dt(VSDTC)) %>%
  filter(!is.na(VSDT)) %>%
  filter(!(is.na(VSSTRESN) & is.na(VSSTRESC))) %>%
  group_by(STUDYID, USUBJID) %>%
  summarise(LSTAVLDT_VS = max(VSDT), .groups = "drop")

# (2) Last complete AE onset date
ae_dates <- ae %>%
  mutate(AEDT = convert_dtc_to_dt(AESTDTC)) %>%
  filter(!is.na(AEDT)) %>%
  group_by(STUDYID, USUBJID) %>%
  summarise(LSTAVLDT_AE = max(AEDT), .groups = "drop")

# (3) Last complete DS date
ds_dates <- ds %>%
  mutate(DSDT = convert_dtc_to_dt(DSSTDTC)) %>%
  filter(!is.na(DSDT)) %>%
  group_by(STUDYID, USUBJID) %>%
  summarise(LSTAVLDT_DS = max(DSDT), .groups = "drop")

# -----------------------------------------------------------------------------
# 10. Derive LSTAVLDT
# -----------------------------------------------------------------------------
adsl <- adsl %>%
  left_join(vs_dates, by = c("STUDYID", "USUBJID")) %>%
  left_join(ae_dates, by = c("STUDYID", "USUBJID")) %>%
  left_join(ds_dates, by = c("STUDYID", "USUBJID")) %>%
  mutate(
    LSTAVLDT_EX = as.Date(TRTEDTM)
  ) %>%
  rowwise() %>%
  mutate(
    .lst_num = {
      vals <- c_across(c(LSTAVLDT_VS, LSTAVLDT_AE, LSTAVLDT_DS, LSTAVLDT_EX))
      vals_num <- as.numeric(vals)
      if (all(is.na(vals_num))) NA_real_ else max(vals_num, na.rm = TRUE)
    },
    LSTAVLDT = as.Date(.lst_num, origin = "1970-01-01")
  ) %>%
  ungroup() %>%
  select(-LSTAVLDT_VS, -LSTAVLDT_AE, -LSTAVLDT_DS, -LSTAVLDT_EX, -.lst_num)

cat("\n=== LSTAVLDT sample ===\n")
print(head(adsl %>% select(USUBJID, TRTSDTM, TRTEDTM, LSTAVLDT), 5))

# -----------------------------------------------------------------------------
# 11. Final ADSL
# -----------------------------------------------------------------------------
adsl_final <- adsl %>%
  arrange(USUBJID) %>%
  select(
    STUDYID, USUBJID, SUBJID, SITEID,
    AGE, AGEU, AGEGR9, AGEGR9N,
    SEX, RACE, ETHNIC, COUNTRY,
    ARM, ACTARM, ARMCD, ACTARMCD,
    RFSTDTC, RFENDTC,
    TRTSDTM, TRTSTMF, TRTEDTM,
    ITTFL,
    LSTAVLDT,
    DTHFL, DTHDTC,
    everything()
  )

# -----------------------------------------------------------------------------
# 12. Validation
# -----------------------------------------------------------------------------
cat("\n=== Final ADSL ===\n")
cat("Dimensions:", nrow(adsl_final), "rows x", ncol(adsl_final), "columns\n")

cat("\n=== AGEGR9 Distribution ===\n")
print(table(adsl_final$AGEGR9, useNA = "always"))

cat("\n=== ITTFL Distribution ===\n")
print(table(adsl_final$ITTFL, useNA = "always"))

cat("\n=== TRTSDTM missingness ===\n")
cat("Non-missing TRTSDTM:", sum(!is.na(adsl_final$TRTSDTM)), "\n")

cat("\n=== LSTAVLDT missingness ===\n")
cat("Non-missing LSTAVLDT:", sum(!is.na(adsl_final$LSTAVLDT)), "\n")

if (nrow(adsl_final) == dplyr::n_distinct(adsl_final$USUBJID)) {
  cat("\n✓ ADSL has exactly one record per subject\n")
} else {
  warning("ADSL has duplicate subjects - review derivations")
}

# -----------------------------------------------------------------------------
# 13. Save outputs
# -----------------------------------------------------------------------------
saveRDS(adsl_final, "adsl.rds")
write.csv(adsl_final, "adsl.csv", row.names = FALSE)

capture.output(
  {
    cat("ADSL creation log\n")
    cat("=================\n")
    cat("Rows:", nrow(adsl_final), "\n")
    cat("Cols:", ncol(adsl_final), "\n\n")
    cat("Variable names:\n")
    print(names(adsl_final))
    cat("\nAGEGR9 distribution:\n")
    print(table(adsl_final$AGEGR9, useNA = "always"))
    cat("\nITTFL distribution:\n")
    print(table(adsl_final$ITTFL, useNA = "always"))
    cat("\nTRTSDTM non-missing:\n")
    print(sum(!is.na(adsl_final$TRTSDTM)))
    cat("\nLSTAVLDT non-missing:\n")
    print(sum(!is.na(adsl_final$LSTAVLDT)))
    cat("\nFirst 10 rows:\n")
    print(head(adsl_final, 10))
  },
  file = "question_2_log.txt"
)

cat("\n✓ ADSL saved to adsl.rds\n")
cat("✓ ADSL saved to adsl.csv\n")
cat("✓ Log saved to question_2_log.txt\n")
cat("\n=== Script completed successfully ===\n")

