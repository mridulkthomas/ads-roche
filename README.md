# Roche PD Data Science – Analytical Data Science Programmer Coding Assessment

> **Candidate Submission** | Roche / Pharmaverse & Python Assessment

---

## Repository Structure

```
.
├── README.md
├── question_1_sdtm/
│   ├── 01_create_ds_domain.R      # SDTM DS domain creation script
│   ├── ds_domain.rds              # Output: DS domain (R native format)
│   ├── ds_domain.csv              # Output: DS domain (CSV for inspection)
│   └── question_1_log.txt         # Output: Script log
│
├── question_2_adam/
│   ├── create_adsl.R              # ADaM ADSL creation script
│   ├── adsl.rds                   # Output: ADSL dataset (R native format)
│   ├── adsl.csv                   # Output: ADSL dataset (CSV for inspection)
│   └── question_2_log.txt         # Output: Script log
│
├── question_3_tlg/
│   ├── 01_create_ae_summary_table.R   # AE summary table (gtsummary)
│   ├── 02_create_visualizations.R     # AE plots (ggplot2)
│   ├── ae_summary_table.html          # Output: Summary table
│   ├── ae_severity_plot.png           # Output: Plot 1 – Severity distribution
│   ├── ae_top10_plot.png              # Output: Plot 2 – Top 10 AEs with CI
│   ├── question_3a_log.txt         # Output: Script log
│   └── question_3b_log.txt         # Output: Script log
│
└── question_4_python/
    └── question_4.py              # GenAI Clinical Data Assistant (Python)
```

---

## Question Summaries

### Question 1 – SDTM DS Domain (`question_1_sdtm/`)

**Objective:** Create an SDTM Disposition (DS) domain from raw clinical trial data using `{sdtm.oak}`.

**Key design decisions:**
- Used `sdtm.oak::assign_no_ct()` for variables without controlled terminology (STUDYID, USUBJID, dates)
- Used `sdtm.oak::assign_ct()` for DSDECOD to map raw collected values to CDISC standard terms via `study_ct`
- DSCAT is derived by context: screening records → `"PROTOCOL MILESTONE"`, others → `"DISPOSITION EVENT"`
- DSSEQ is generated last, after sorting by USUBJID and DSSTDTC, to ensure stable sequence numbers
- DSSTDY is calculated using the standard SDTM study-day formula: `DSSTDT − RFSTDT + 1` (positive days), `DSSTDT − RFSTDT` (negative days)

**Output variables:** `STUDYID, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD, DSCAT, VISITNUM, VISIT, DSDTC, DSSTDTC, DSSTDY`

**R packages:** `sdtm.oak`, `pharmaverseraw`, `dplyr`, `lubridate`

---

### Question 2 – ADaM ADSL Dataset (`question_2_adam/`)

**Objective:** Create a Subject-Level Analysis Dataset (ADSL) using `{admiral}`.

**Key design decisions:**

| Variable | Derivation Approach |
|---|---|
| `AGEGR9` / `AGEGR9N` | `dplyr::case_when()` with age thresholds `<18`, `18–50`, `>50` mapped to numeric codes 1, 2, 3 |
| `TRTSDTM` / `TRTSTMF` | `admiral::derive_vars_dtm()` on `EX.EXSTDTC` with `highest_imputation = "M"` (imputes hours/minutes but not seconds) |
| `ITTFL` | `"Y"` when `DM.ARM` is not missing; `"N"` otherwise |
| `LSTAVLDT` | `rowwise()` maximum across 4 source dates: last VS date (with valid result), last AE start date, last DS date, and datepart of `TRTEDTM` |

**Valid dose filter (for TRTSDTM/TRTEDTM):** `EXDOSE > 0` OR (`EXDOSE == 0` AND `EXTRT` contains `"PLACEBO"`)

**R packages:** `admiral`, `pharmaversesdtm`, `dplyr`, `lubridate`, `stringr`, `tidyr`

---

### Question 3 – TLG: Adverse Events Reporting (`question_3_tlg/`)

**Objective:** Create regulatory-style AE tables and visualizations.

#### Task 1 – Summary Table (`01_create_ae_summary_table.R`)
- Filters to TEAEs (`TRTEMFL == "Y"`)
- Uses `{gtsummary}` `tbl_summary()` stratified by `ACTARM`
- Adds an overall/total column via `add_overall()`
- Denominators are the full Safety Population from ADSL (not just subjects with AEs)
- Output: `ae_summary_table.html`

#### Task 2 – Visualizations (`02_create_visualizations.R`)

**Plot 1 – AE Severity Distribution (stacked bar chart):**
- Counts AE events (not subjects) by `ACTARM` and `AESEV`
- Ordered factor: `MILD → MODERATE → SEVERE`
- Custom colour palette matching sample output
- Output: `ae_severity_plot.png`

**Plot 2 – Top 10 Most Frequent AEs with 95% CI:**
- Subject-level incidence (counts distinct `USUBJID` per `AETERM`, not events)
- Clopper-Pearson exact binomial 95% CIs via `{binom}` package
- Horizontal dot-and-whisker plot ordered by frequency
- Output: `ae_top10_plot.png`

**R packages:** `ggplot2`, `gtsummary`, `gt`, `pharmaverseadam`, `dplyr`, `forcats`, `binom`, `scales`

---

### Question 4 – GenAI Clinical Data Assistant (`question_4_python/`)

**Objective:** Build an LLM-powered agent that translates natural language questions into Pandas queries.

**Architecture – `ClinicalTrialDataAgent`:**

```
User Question
     │
     ▼
┌─────────────────────────────────────┐
│  parse_question()                   │
│  ┌──────────────────────────────┐   │
│  │  System Prompt               │   │
│  │  + ADAE Schema (JSON)        │   │──► LLM (OpenAI / Anthropic / Mock)
│  │  + User Question             │   │
│  └──────────────────────────────┘   │
│  → Returns: {target_column,         │
│              filter_value}          │
└─────────────────────────────────────┘
     │
     ▼
┌─────────────────────────────────────┐
│  execute_filter()                   │
│  Applies Pandas boolean filter on   │
│  the ADAE DataFrame                 │
└─────────────────────────────────────┘
     │
     ▼
  QueryResult
  ├── n_subjects (unique USUBJID count)
  ├── subject_ids (list)
  └── n_records
```

**Key design decisions:**
- **Schema-driven prompting:** A full JSON schema of ADAE columns (including `maps_to` synonyms) is injected into every system prompt, enabling the LLM to handle synonymous phrasing (e.g., "intensity" → `AESEV`)
- **Structured JSON output:** The prompt instructs the LLM to return only a JSON object with `target_column` and `filter_value` – no prose, no markdown
- **Graceful degradation:** Three LLM backends in priority order: OpenAI (via LangChain) → Anthropic Claude → mock/rule-based fallback
- **Case-insensitive matching:** Pandas filter uses `.str.upper()` comparison to handle casing inconsistencies in raw data
- **JSON extraction with regex:** `re.search(r"\{.*?\}", ...)` handles cases where an LLM wraps JSON in markdown code fences

**Python dependencies:** `pandas`, `langchain-openai` (optional), `anthropic` (optional)

**Usage:**
```python
from question_4 import ClinicalTrialDataAgent

# With OpenAI API key
agent = ClinicalTrialDataAgent(api_key="sk-...")

# Without API key (mock mode – demonstrates full Prompt → Parse → Execute flow)
agent = ClinicalTrialDataAgent()

result = agent.query("Give me subjects with moderate severity AEs")
print(result)
# → n_subjects: 47, subject_ids: [...], filter: {target_column: AESEV, filter_value: MODERATE}
```

---

## Setup & Running

### R (Questions 1–3)

```r
# Install required packages
install.packages(c(
  "admiral", "sdtm.oak", "pharmaverseraw", "pharmaversesdtm", "pharmaverseadam",
  "gt", "gtsummary", "ggplot2", "dplyr", "tidyr", "lubridate",
  "stringr", "forcats", "scales", "binom"
))

# Run scripts from the repo root directory
source("question_1_sdtm/01_create_ds_domain.R")
source("question_2_adam/create_adsl.R")
source("question_3_tlg/01_create_ae_summary_table.R")
source("question_3_tlg/02_create_visualizations.R")
```

**R version:** 4.2.0 or above  
**Environment:** Works on local R installation or Posit Cloud (free tier sufficient)

### Python (Question 4)

```bash
# Install dependencies
pip install pandas

# Optional: for real LLM support
pip install langchain-openai     # OpenAI via LangChain
pip install anthropic             # Anthropic Claude

# Run test script
python question_4_python/question_4.py
```


---
