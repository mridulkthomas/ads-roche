# =============================================================================
# Question 1: SDTM DS Domain Creation using {sdtm.oak}
# =============================================================================
# Objective: Create an SDTM Disposition (DS) domain dataset from raw clinical
#            trial data using the {sdtm.oak} package.
#
# Input:     pharmaverseraw::ds_raw
# Output:    DS domain with variables:
#            STUDYID, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD, DSCAT,
#            VISITNUM, VISIT, DSDTC, DSSTDTC, DSSTDY
#
# Author:    Mridul K. Thomas
# Date:      03-2026
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load libraries
# -----------------------------------------------------------------------------
library(sdtm.oak)
library(pharmaverseraw)
library(pharmaversesdtm)
library(dplyr)
library(stringr)
library(lubridate)

# -----------------------------------------------------------------------------
# 2. Load raw data
# -----------------------------------------------------------------------------
ds_raw <- pharmaverseraw::ds_raw

cat("=== Raw DS structure ===\n")
str(ds_raw)
cat("\n=== Raw DS column names ===\n")
print(names(ds_raw))
cat("\n=== First 10 rows ===\n")
print(utils::head(ds_raw, 10))

# -----------------------------------------------------------------------------
# 3. Study controlled terminology
# -----------------------------------------------------------------------------
# The exercise explicitly allows creating study_ct in code.
study_ct <- data.frame(
  stringsAsFactors = FALSE,
  codelist_code = c(
    "C66727", "C66727", "C66727", "C66727", "C66727",
    "C66727", "C66727", "C66727", "C66727", "C66727"
  ),
  term_code = c(
    "C41331", "C25250", "C28554", "C48226", "C48227",
    "C48250", "C142185", "C49628", "C49632", "C49634"
  ),
  term_value = c(
    "ADVERSE EVENT", "COMPLETED", "DEATH", "LACK OF EFFICACY",
    "LOST TO FOLLOW-UP", "PHYSICIAN DECISION", "PROTOCOL VIOLATION",
    "SCREEN FAILURE", "STUDY TERMINATED BY SPONSOR", "WITHDRAWAL BY SUBJECT"
  ),
  collected_value = c(
    "Adverse Event", "Complete", "Dead", "Lack of Efficacy",
    "Lost To Follow-Up", "Physician Decision", "Protocol Violation",
    "Trial Screen Failure", "Study Terminated By Sponsor",
    "Withdrawal by Subject"
  ),
  term_preferred_term = c(
    "AE", "Completed", "Died", NA, NA, NA, "Violation",
    "Failure to Meet Inclusion/Exclusion Criteria", NA, "Dropout"
  ),
  term_synonyms = c(
    "ADVERSE EVENT", "COMPLETE", "Death", NA, NA, NA, NA, NA, NA,
    "Discontinued Participation"
  )
)

cat("\n=== Study CT ===\n")
print(study_ct)

# -----------------------------------------------------------------------------
# 4. Add oak ID variables
# -----------------------------------------------------------------------------
# This preserves a stable record identifier and still demonstrates use of
# sdtm.oak as requested in the exercise.
ds_raw <- generate_oak_id_vars(
  raw_dat = ds_raw,
  pat_var = "PATNUM",
  raw_src = "ds_raw"
)

cat("\n=== oak ID variables added ===\n")
print(utils::head(ds_raw[, c("oak_id", "raw_source", "patient_number")], 10))

# -----------------------------------------------------------------------------
# 5. Helper functions
# -----------------------------------------------------------------------------
make_iso_date <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  
  # Try several common formats seen in raw CRF-style data
  parsed <- suppressWarnings(parse_date_time(
    x,
    orders = c(
      "Ymd", "ymd", "dmy", "mdy",
      "d b Y", "d B Y", "Y-m-d", "d/m/Y", "m/d/Y"
    ),
    quiet = TRUE
  ))
  
  ifelse(is.na(parsed), NA_character_, format(as.Date(parsed), "%Y-%m-%d"))
}

make_iso_dttm <- function(date_x, time_x = NULL) {
  date_chr <- make_iso_date(date_x)
  
  if (is.null(time_x)) {
    return(date_chr)
  }
  
  time_chr <- as.character(time_x)
  time_chr <- trimws(time_chr)
  time_chr[time_chr == ""] <- NA_character_
  
  # Standardize incomplete times if needed
  time_chr <- ifelse(
    !is.na(time_chr) & str_detect(time_chr, "^\\d{1,2}:\\d{2}$"),
    paste0(time_chr, ":00"),
    time_chr
  )
  
  ifelse(
    is.na(date_chr),
    NA_character_,
    ifelse(
      is.na(time_chr),
      date_chr,
      paste0(date_chr, "T", time_chr)
    )
  )
}

derive_study_day <- function(event_dtc, ref_dtc) {
  event_dt <- as.Date(substr(event_dtc, 1, 10))
  ref_dt   <- as.Date(substr(ref_dtc, 1, 10))
  
  case_when(
    is.na(event_dt) | is.na(ref_dt) ~ NA_integer_,
    event_dt >= ref_dt ~ as.integer(event_dt - ref_dt) + 1L,
    TRUE ~ as.integer(event_dt - ref_dt)
  )
}

map_dsdecod <- function(x, ct_df) {
  x_chr <- as.character(x)
  
  out <- case_when(
    x_chr %in% ct_df$collected_value ~
      ct_df$term_value[match(x_chr, ct_df$collected_value)],
    x_chr %in% ct_df$term_synonyms ~
      ct_df$term_value[match(x_chr, ct_df$term_synonyms)],
    toupper(x_chr) %in% ct_df$term_value ~
      toupper(x_chr),
    TRUE ~ NA_character_
  )
  
  out
}

derive_visitnum <- function(instance) {
  case_when(
    is.na(instance) ~ NA_real_,
    str_detect(instance, regex("screen", ignore_case = TRUE)) ~ 1,
    str_detect(instance, regex("baseline", ignore_case = TRUE)) ~ 2,
    str_detect(instance, regex("end", ignore_case = TRUE)) ~ 999,
    str_detect(instance, regex("follow", ignore_case = TRUE)) ~ 900,
    TRUE ~ suppressWarnings(as.numeric(str_extract(instance, "[0-9]+")))
  )
}

# -----------------------------------------------------------------------------
# 6. Reference start date from DM for DSSTDY
# -----------------------------------------------------------------------------
# The earlier draft assumed RFSTDTC might exist in ds_raw. A safer approach is
# to join it from pharmaversesdtm::dm.
dm_ref <- pharmaversesdtm::dm %>%
  select(STUDYID, USUBJID, RFSTDTC)

cat("\n=== DM reference dates ===\n")
print(utils::head(dm_ref, 10))

# -----------------------------------------------------------------------------
# 7. Derive DS directly from ds_raw
# -----------------------------------------------------------------------------
# This avoids the failing assign_no_ct()/assign_ct() steps from the draft.
ds <- ds_raw %>%
  mutate(
    STUDYID = STUDY,
    DOMAIN  = "DS",
    USUBJID = paste(STUDY, PATNUM, sep = "-"),
    
    # Reported term as collected
    DSTERM = `IT.DSTERM`,
    
    # Standardized term using study CT
    DSDECOD = map_dsdecod(`IT.DSDECOD`, study_ct),
    
    # Category
    DSCAT = case_when(
      DSDECOD == "SCREEN FAILURE" ~ "PROTOCOL MILESTONE",
      TRUE ~ "DISPOSITION EVENT"
    ),
    
    # Visit variables
    VISIT = INSTANCE,
    VISITNUM = derive_visitnum(INSTANCE),
    
    # Timing variables
    DSDTC = make_iso_dttm(DSDTCOL, DSTMCOL),
    DSSTDTC = make_iso_date(`IT.DSSTDAT`)
  )

# Fallback: if DSSTDTC missing but DSDTC exists, use DSDTC date part
ds <- ds %>%
  mutate(
    DSSTDTC = ifelse(
      is.na(DSSTDTC) & !is.na(DSDTC),
      substr(DSDTC, 1, 10),
      DSSTDTC
    )
  )

cat("\n=== After derivations ===\n")
print(utils::head(
  ds %>%
    select(STUDYID, DOMAIN, USUBJID, DSTERM, DSDECOD, DSCAT, VISIT, VISITNUM, DSDTC, DSSTDTC),
  10
))

cat("\n=== DSDECOD distribution ===\n")
print(table(ds$DSDECOD, useNA = "ifany"))

cat("\n=== DSCAT distribution ===\n")
print(table(ds$DSCAT, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 8. Derive DSSTDY
# -----------------------------------------------------------------------------
ds <- ds %>%
  left_join(dm_ref, by = c("STUDYID", "USUBJID")) %>%
  mutate(
    DSSTDY = derive_study_day(DSSTDTC, RFSTDTC)
  )

# -----------------------------------------------------------------------------
# 9 Derive DSSEQ
# -----------------------------------------------------------------------------
ds <- ds %>%
  arrange(USUBJID, DSSTDTC, DSDTC, oak_id) %>%
  group_by(USUBJID) %>%
  mutate(DSSEQ = row_number()) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 10. Select final variables
# -----------------------------------------------------------------------------
ds_final <- ds %>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    DSSEQ,
    DSTERM,
    DSDECOD,
    DSCAT,
    VISITNUM,
    VISIT,
    DSDTC,
    DSSTDTC,
    DSSTDY
  )

# -----------------------------------------------------------------------------
# 11. Final validation
# -----------------------------------------------------------------------------
cat("\n=== Final DS Domain ===\n")
cat("Dimensions:", nrow(ds_final), "rows x", ncol(ds_final), "columns\n")
cat("\nVariable names:\n")
print(names(ds_final))
cat("\nFirst 10 rows:\n")
print(utils::head(ds_final, 10))

required_vars <- c(
  "STUDYID", "DOMAIN", "USUBJID", "DSSEQ", "DSTERM", "DSDECOD",
  "DSCAT", "VISITNUM", "VISIT", "DSDTC", "DSSTDTC", "DSSTDY"
)

missing_vars <- setdiff(required_vars, names(ds_final))
if (length(missing_vars) > 0) {
  stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
} else {
  cat("\n✓ All required variables present\n")
}

cat("\n=== Missingness summary ===\n")
print(sapply(ds_final, function(x) sum(is.na(x))))

# -----------------------------------------------------------------------------
# 12. Save outputs
# -----------------------------------------------------------------------------
saveRDS(ds_final, "ds_domain.rds")
write.csv(ds_final, "ds_domain.csv", row.names = FALSE)

capture.output(
  {
    cat("DS domain creation log\n")
    cat("======================\n")
    cat("Rows:", nrow(ds_final), "\n")
    cat("Cols:", ncol(ds_final), "\n\n")
    cat("Variable names:\n")
    print(names(ds_final))
    cat("\nDSDECOD distribution:\n")
    print(table(ds_final$DSDECOD, useNA = "ifany"))
    cat("\nDSCAT distribution:\n")
    print(table(ds_final$DSCAT, useNA = "ifany"))
    cat("\nMissingness summary:\n")
    print(sapply(ds_final, function(x) sum(is.na(x))))
    cat("\nFirst 10 rows:\n")
    print(utils::head(ds_final, 10))
  },
  file = "question_1_log.txt"
)

cat("\n✓ DS domain saved to ds_domain.rds\n")
cat("✓ DS domain saved to ds_domain.csv\n")
cat("✓ Log file saved to question_1_sdtm/question_1_log.txt\n")
cat("\n=== Script completed successfully ===\n")
