# =============================================================================
# Question 3b: TLG - Adverse Events Visualizations using {ggplot2}
# =============================================================================
# Objective: Create two ggplot2 visualizations for adverse events reporting.
#
# Plot 1: AE severity distribution by treatment arm (stacked bar chart)
#         - Variable: AESEV (MILD, MODERATE, SEVERE)
#         - X-axis: Treatment Arm
#         - Y-axis: Count of AEs
#
# Plot 2: Top 10 most frequent AEs with 95% CI for incidence rates
#         - Variable: AETERM
#         - X-axis: Percentage of Patients (%)
#         - CI method: Clopper-Pearson (exact binomial)
#
# Input:  pharmaverseadam::adae, pharmaverseadam::adsl
# Output: ae_severity_plot.png, ae_top10_plot.png
#
# Author:    Mridul K. Thomas
# Date:      03-2026
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load Required Libraries
# -----------------------------------------------------------------------------
library(ggplot2)         # Visualization
library(pharmaverseadam) # ADaM example datasets
library(dplyr)           # Data manipulation
library(tidyr)           # Data reshaping
library(forcats)         # Factor reordering
library(scales)          # Axis formatting
library(binom)           # Exact binomial CI (Clopper-Pearson)

# -----------------------------------------------------------------------------
# 2. Load Data
# -----------------------------------------------------------------------------
adae <- pharmaverseadam::adae
adsl <- pharmaverseadam::adsl

# Filter to Treatment-Emergent AEs only
adae_te <- adae %>%
  filter(TRTEMFL == "Y")

# Get total N per treatment arm for denominators
n_by_arm <- adsl %>%
  count(ACTARM, name = "N")

n_total <- nrow(adsl)

cat("=== Data loaded ===\n")
cat("TEAE records:", nrow(adae_te), "\n")
cat("Total subjects:", n_total, "\n")

# =============================================================================
# PLOT 1: AE Severity Distribution by Treatment Arm
# =============================================================================

# -----------------------------------------------------------------------------
# 3. Prepare Data for Plot 1
# -----------------------------------------------------------------------------
# Count AEs by treatment arm and severity
severity_data <- adae_te %>%
  filter(!is.na(AESEV)) %>%
  count(ACTARM, AESEV, name = "count") %>%
  # Ensure severity is an ordered factor for consistent stacking
  mutate(
    AESEV = factor(AESEV, levels = c("MILD", "MODERATE", "SEVERE"))
  ) %>%
  # Add arm N for labeling
  left_join(n_by_arm, by = "ACTARM") %>%
  mutate(
    arm_label = sprintf("%s\n(N = %d)", ACTARM, N)
  )

cat("\n=== Severity Data ===\n")
print(severity_data)

# Define colour palette matching sample output
severity_colours <- c(
  "MILD"     = "#E8735A",   # Salmon/orange-red
  "MODERATE" = "#2E8B57",   # Sea green
  "SEVERE"   = "#6495ED"    # Cornflower blue
)

# -----------------------------------------------------------------------------
# 4. Create Plot 1
# -----------------------------------------------------------------------------
plot1 <- ggplot(
  severity_data,
  aes(x = arm_label, y = count, fill = AESEV)
) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(
    values = severity_colours,
    name   = "Severity/Intensity",
    breaks = c("MILD", "MODERATE", "SEVERE")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = seq(0, 500, by = 100)
  ) +
  labs(
    title   = "AE severity distribution by treatment",
    x       = "Treatment Arm",
    y       = "Count of AEs",
    caption = "Treatment-emergent adverse events only (TRTEMFL = 'Y')"
  ) +
  theme_classic(base_size = 12) +
  theme(
    # Title
    plot.title          = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.caption        = element_text(size = 9, colour = "grey50"),
    # Axes
    axis.title          = element_text(face = "bold"),
    axis.text.x         = element_text(size = 10),
    axis.text.y         = element_text(size = 10),
    # Legend
    legend.title        = element_text(face = "bold"),
    legend.position     = "right",
    legend.key.size     = unit(0.5, "cm"),
    # Grid
    panel.grid.major.y  = element_line(colour = "grey90", linetype = "dashed"),
    panel.grid.major.x  = element_blank(),
    # Background
    plot.background     = element_rect(fill = "white", colour = NA),
    panel.background    = element_rect(fill = "white", colour = NA)
  )

# Save Plot 1
ggsave(
  filename = "ae_severity_plot.png",
  plot     = plot1,
  width    = 8,
  height   = 6,
  dpi      = 300,
  bg       = "white"
)

cat("\n✓ Plot 1 (AE Severity) saved to ae_severity_plot.png\n")
print(plot1)

# =============================================================================
# PLOT 2: Top 10 Most Frequent AEs with 95% Clopper-Pearson CI
# =============================================================================

# -----------------------------------------------------------------------------
# 5. Prepare Data for Plot 2
# -----------------------------------------------------------------------------
# Count unique subjects with each AE (not AE events - subject-level incidence)
ae_incidence <- adae_te %>%
  distinct(USUBJID, AETERM) %>%        # One record per subject per AE term
  count(AETERM, name = "n_subjects") %>%
  arrange(desc(n_subjects)) %>%
  slice_head(n = 10) %>%               # Top 10 most frequent
  mutate(
    N       = n_total,
    # Point estimate: incidence rate as percentage
    pct     = n_subjects / N * 100,
    # Clopper-Pearson exact 95% CI
    ci_lo   = binom.confint(n_subjects, N, methods = "exact")$lower * 100,
    ci_hi   = binom.confint(n_subjects, N, methods = "exact")$upper * 100,
    # Reorder factor by frequency (highest at top in horizontal plot)
    AETERM  = fct_reorder(AETERM, pct)
  )

cat("\n=== Top 10 AEs ===\n")
print(ae_incidence %>% select(AETERM, n_subjects, pct, ci_lo, ci_hi))

# -----------------------------------------------------------------------------
# 6. Create Plot 2
# -----------------------------------------------------------------------------
plot2 <- ggplot(
  ae_incidence,
  aes(x = pct, y = AETERM)
) +
  # Confidence interval lines
  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi),
    height  = 0.3,
    size    = 0.8,
    colour  = "black"
  ) +
  # Point estimate
  geom_point(
    size   = 3,
    colour = "black",
    fill   = "black"
  ) +
  # Vertical reference line at 0
  geom_vline(xintercept = 0, colour = "grey70", linetype = "solid") +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, NA),
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  labs(
    title    = "Top 10 Most Frequent Adverse Events",
    subtitle = sprintf("n = %d subjects; 95%% Clopper-Pearson CIs", n_total),
    x        = "Percentage of Patients (%)",
    y        = NULL,
    caption  = "Treatment-emergent adverse events only (TRTEMFL = 'Y'); Subject-level incidence"
  ) +
  theme_classic(base_size = 12) +
  theme(
    # Title
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, colour = "grey40"),
    plot.caption     = element_text(size = 9, colour = "grey50"),
    # Axes
    axis.title.x     = element_text(face = "bold"),
    axis.text.y      = element_text(size = 10, hjust = 1),
    axis.text.x      = element_text(size = 10),
    # Grid
    panel.grid.major.x = element_line(colour = "grey90", linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    # Background
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# Save Plot 2
ggsave(
  filename = "ae_top10_plot.png",
  plot     = plot2,
  width    = 9,
  height   = 6,
  dpi      = 300,
  bg       = "white"
)

cat("\n✓ Plot 2 (Top 10 AEs) saved to ae_top10_plot.png\n")
print(plot2)

capture.output(
  {
    cat("Question 3b log\n")
    cat("================\n")
    cat("ADAE rows:", nrow(adae), "\n")
    cat("ADSL rows:", nrow(adsl), "\n")
    cat("TEAE rows:", nrow(adae_te), "\n")
    cat("Total subjects:", n_total, "\n\n")
    
    cat("Subjects per arm:\n")
    print(n_by_arm)
    
    cat("\nSeverity data:\n")
    print(severity_data)
    
    cat("\nTop 10 AE incidence data:\n")
    print(ae_incidence %>% select(AETERM, n_subjects, pct, ci_lo, ci_hi))
    
    cat("\nOutput files:\n")
    cat("ae_severity_plot.png\n")
    cat("ae_top10_plot.png\n")
  },
  file = "question_3b_log.txt"
)

cat("\n✓ Log saved to question_3b_log.txt\n")
cat("\n=== Script completed successfully ===\n")
