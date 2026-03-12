"""
=============================================================================
Question 4: GenAI Clinical Data Assistant (LLM & LangChain)
=============================================================================
Objective: A Generative AI assistant that translates natural language questions
           into structured Pandas queries against a clinical AE dataset.

The agent:
  1. Understands the ADAE dataset schema
  2. Uses an LLM to parse user questions into JSON (target_column, filter_value)
  3. Executes the resulting Pandas filter
  4. Returns count of unique subjects (USUBJID) and their IDs

Author:    Candidate
Date:      2024

Usage:
  # With real API key:
  agent = ClinicalTrialDataAgent(api_key="sk-...")

  # Without API key (mock mode):
  agent = ClinicalTrialDataAgent(api_key=None)

  result = agent.query("Show me subjects with moderate severity AEs")
  print(result)
=============================================================================
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from typing import Any

import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Optional imports – gracefully degrade when packages are not installed
# ─────────────────────────────────────────────────────────────────────────────
try:
    from langchain_openai import ChatOpenAI
    from langchain_core.prompts import ChatPromptTemplate
    from langchain_core.output_parsers import StrOutputParser
    LANGCHAIN_AVAILABLE = True
except ImportError:
    LANGCHAIN_AVAILABLE = False
    print("[INFO] LangChain / OpenAI packages not found. Running in mock mode.")

try:
    import anthropic as anthropic_sdk
    ANTHROPIC_AVAILABLE = True
except ImportError:
    ANTHROPIC_AVAILABLE = False


# =============================================================================
# 1. Dataset Schema Definition
# =============================================================================
# This schema description is injected into every LLM prompt so the model
# knows exactly which columns exist and what values they contain.

ADAE_SCHEMA: dict[str, dict[str, Any]] = {
    "STUDYID": {
        "description": "Unique study identifier",
        "type": "string",
        "example": "CDISCPILOT01",
    },
    "USUBJID": {
        "description": "Unique subject identifier (primary key per subject)",
        "type": "string",
        "example": "CDISCPILOT01-01-701-1015",
    },
    "AESOC": {
        "description": (
            "Primary System Organ Class (SOC) – high-level body system classification. "
            "Examples: 'CARDIAC DISORDERS', 'SKIN AND SUBCUTANEOUS TISSUE DISORDERS', "
            "'NERVOUS SYSTEM DISORDERS', 'GASTROINTESTINAL DISORDERS'"
        ),
        "type": "string",
        "maps_to": ["body system", "organ class", "system organ class", "SOC"],
    },
    "AETERM": {
        "description": (
            "Reported (verbatim) adverse event term – specific condition name. "
            "Examples: 'HEADACHE', 'NAUSEA', 'DIZZINESS', 'APPLICATION SITE PRURITUS'"
        ),
        "type": "string",
        "maps_to": ["adverse event", "condition", "specific AE", "term", "diagnosis"],
    },
    "AEDECOD": {
        "description": (
            "Dictionary-derived AE term (MedDRA preferred term). "
            "Similar to AETERM but standardized."
        ),
        "type": "string",
        "maps_to": ["preferred term", "MedDRA term"],
    },
    "AESEV": {
        "description": (
            "AE severity / intensity grade. "
            "Allowed values: 'MILD', 'MODERATE', 'SEVERE'"
        ),
        "type": "categorical",
        "allowed_values": ["MILD", "MODERATE", "SEVERE"],
        "maps_to": ["severity", "intensity", "grade", "how bad", "serious"],
    },
    "AESER": {
        "description": "Serious AE flag. Allowed values: 'Y' (serious) or 'N' (not serious).",
        "type": "categorical",
        "allowed_values": ["Y", "N"],
        "maps_to": ["serious", "SAE", "serious adverse event"],
    },
    "ACTARM": {
        "description": (
            "Actual treatment arm. "
            "Examples: 'Placebo', 'Xanomeline High Dose', 'Xanomeline Low Dose'"
        ),
        "type": "string",
        "maps_to": ["treatment", "arm", "dose group", "treatment group"],
    },
    "TRTEMFL": {
        "description": "Treatment-emergent AE flag. Allowed values: 'Y' or 'N'.",
        "type": "categorical",
        "allowed_values": ["Y", "N"],
        "maps_to": ["treatment-emergent", "TEAE", "emergent"],
    },
    "AESTDTC": {
        "description": "AE start date in ISO 8601 format (YYYY-MM-DD).",
        "type": "date",
        "maps_to": ["start date", "onset date", "when did it start"],
    },
    "AEENDTC": {
        "description": "AE end date in ISO 8601 format (YYYY-MM-DD).",
        "type": "date",
        "maps_to": ["end date", "resolution date", "when did it end"],
    },
}

SCHEMA_STR = json.dumps(
    {k: {"description": v["description"], "type": v.get("type")}
     for k, v in ADAE_SCHEMA.items()},
    indent=2,
)


# =============================================================================
# 2. Prompt Template
# =============================================================================

SYSTEM_PROMPT = """You are a clinical data analyst assistant with expertise in CDISC ADaM datasets.
Your task is to translate a natural language question about adverse events into a structured JSON filter.

Here is the schema of the ADAE (Adverse Events Analysis Dataset):
{schema}

IMPORTANT INSTRUCTIONS:
1. Respond ONLY with a valid JSON object. Do not include any explanation, markdown, or extra text.
2. The JSON must have exactly these two keys:
   - "target_column": The ADAE column name to filter on (must be one of the schema columns above)
   - "filter_value": The exact value to filter for (use UPPERCASE for categorical variables like AESEV, AESER, TRTEMFL)
3. For AESEV: map "moderate" → "MODERATE", "mild" → "MILD", "severe" → "SEVERE"
4. For AESER: map "serious" → "Y", "not serious" → "N"
5. For TRTEMFL: map "treatment-emergent" → "Y"
6. For AESOC and AETERM: use UPPERCASE as these are stored in uppercase in the dataset
7. If the question is ambiguous, use your best judgment to select the most appropriate column.

Example input: "Show me subjects with moderate adverse events"
Example output: {{"target_column": "AESEV", "filter_value": "MODERATE"}}

Example input: "Which patients had headache?"
Example output: {{"target_column": "AETERM", "filter_value": "HEADACHE"}}

Example input: "Find AEs in the cardiac system"
Example output: {{"target_column": "AESOC", "filter_value": "CARDIAC DISORDERS"}}
"""

USER_PROMPT = "Question: {question}"


# =============================================================================
# 3. Mock LLM (for when no API key is available)
# =============================================================================

def _mock_llm_parse(question: str) -> dict[str, str]:
    """
    Rule-based fallback that simulates LLM output.
    Demonstrates the Prompt -> Parse -> Execute flow without an API key.
    """
    q = question.lower()

    # Severity mapping
    for sev in ("moderate", "mild", "severe"):
        if sev in q:
            return {"target_column": "AESEV", "filter_value": sev.upper()}

    # Serious AE
    if "serious" in q and "not serious" not in q:
        return {"target_column": "AESER", "filter_value": "Y"}
    if "not serious" in q:
        return {"target_column": "AESER", "filter_value": "N"}

    # Treatment-emergent
    if "treatment-emergent" in q or "emergent" in q or "teae" in q:
        return {"target_column": "TRTEMFL", "filter_value": "Y"}

    # Treatment arm
    for arm_kw, arm_val in [
        ("placebo", "Placebo"),
        ("high dose", "Xanomeline High Dose"),
        ("low dose", "Xanomeline Low Dose"),
    ]:
        if arm_kw in q:
            return {"target_column": "ACTARM", "filter_value": arm_val}

    # Body system / SOC keywords
    soc_map = {
        "cardiac":         "CARDIAC DISORDERS",
        "skin":            "SKIN AND SUBCUTANEOUS TISSUE DISORDERS",
        "nervous":         "NERVOUS SYSTEM DISORDERS",
        "gastro":          "GASTROINTESTINAL DISORDERS",
        "general disorder": "GENERAL DISORDERS AND ADMINISTRATION SITE CONDITIONS",
    }
    for kw, soc in soc_map.items():
        if kw in q:
            return {"target_column": "AESOC", "filter_value": soc}

    # Specific AE terms (common ones)
    ae_terms = [
        "headache", "nausea", "dizziness", "pruritus", "rash",
        "vomiting", "fatigue", "diarrhoea", "erythema", "application site",
    ]
    for term in ae_terms:
        if term in q:
            return {"target_column": "AETERM", "filter_value": term.upper()}

    # Default fallback
    return {"target_column": "TRTEMFL", "filter_value": "Y"}


# =============================================================================
# 4. ClinicalTrialDataAgent Class
# =============================================================================

@dataclass
class QueryResult:
    """Structured result from a clinical data query."""
    question:        str
    parsed_filter:   dict[str, str]
    n_subjects:      int
    subject_ids:     list[str]
    n_records:       int
    llm_mode:        str   # "openai", "anthropic", or "mock"

    def __str__(self) -> str:
        ids_preview = (
            self.subject_ids[:5] + ["..."]
            if len(self.subject_ids) > 5
            else self.subject_ids
        )
        return (
            f"\n{'='*60}\n"
            f"Question    : {self.question}\n"
            f"Filter      : {self.parsed_filter}\n"
            f"LLM Mode    : {self.llm_mode}\n"
            f"N Subjects  : {self.n_subjects}\n"
            f"N Records   : {self.n_records}\n"
            f"Subject IDs : {ids_preview}\n"
            f"{'='*60}"
        )


class ClinicalTrialDataAgent:
    """
    An AI-powered agent that answers natural language questions about
    adverse events by translating them into Pandas DataFrame filters.

    Supports three LLM backends (in priority order):
      1. OpenAI via LangChain  (requires OPENAI_API_KEY or api_key param)
      2. Anthropic Claude      (requires ANTHROPIC_API_KEY)
      3. Mock / rule-based     (no API key needed)

    Parameters
    ----------
    df : pd.DataFrame, optional
        The ADAE dataframe to query. If None, loads from adae.csv.
    api_key : str, optional
        OpenAI API key. Falls back to env var OPENAI_API_KEY.
    anthropic_api_key : str, optional
        Anthropic API key. Falls back to env var ANTHROPIC_API_KEY.
    model : str
        OpenAI model name (default: "gpt-4o-mini").
    """

    def __init__(
        self,
        df: pd.DataFrame | None = None,
        api_key: str | None = None,
        anthropic_api_key: str | None = None,
        model: str = "gpt-4o-mini",
    ) -> None:
        # ── Load data ───────────────────────────────────────────────────────
        if df is not None:
            self.df = df.copy()
        else:
            self.df = self._load_default_data()

        # ── LLM setup ───────────────────────────────────────────────────────
        self._llm_mode = "mock"
        self._chain = None

        openai_key = api_key or os.getenv("OPENAI_API_KEY")
        anth_key   = anthropic_api_key or os.getenv("ANTHROPIC_API_KEY")

        if LANGCHAIN_AVAILABLE and openai_key:
            self._setup_openai_chain(openai_key, model)
            self._llm_mode = "openai"
        elif ANTHROPIC_AVAILABLE and anth_key:
            self._anthropic_key = anth_key
            self._llm_mode = "anthropic"
        else:
            print("[INFO] No API key provided or packages missing. Using mock LLM.")

        print(f"[INFO] ClinicalTrialDataAgent initialized (mode={self._llm_mode})")
        print(f"[INFO] Dataset shape: {self.df.shape}")

    # ──────────────────────────────────────────────────────────────────────
    # Private helpers
    # ──────────────────────────────────────────────────────────────────────

    @staticmethod
    def _load_default_data() -> pd.DataFrame:
        """Load ADAE data from CSV or generate synthetic fallback."""
        csv_path = "adae.csv"
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            print(f"[INFO] Loaded {len(df)} rows from {csv_path}")
            return df

        # Synthetic fallback for demonstration
        print("[INFO] adae.csv not found – generating synthetic ADAE data.")
        import numpy as np
        rng = np.random.default_rng(42)
        n = 500
        subjects  = [f"CDISCPILOT01-01-{rng.integers(700, 715)}-{rng.integers(1000, 1999)}"
                     for _ in range(n)]
        arms      = rng.choice(["Placebo", "Xanomeline High Dose", "Xanomeline Low Dose"], n)
        severities = rng.choice(["MILD", "MODERATE", "SEVERE"], n, p=[0.6, 0.3, 0.1])
        ae_terms  = rng.choice(
            ["HEADACHE", "NAUSEA", "DIZZINESS", "APPLICATION SITE PRURITUS",
             "RASH", "FATIGUE", "VOMITING", "ERYTHEMA", "DIARRHOEA", "PRURITUS"],
            n,
        )
        soc_map = {
            "HEADACHE": "NERVOUS SYSTEM DISORDERS",
            "NAUSEA":   "GASTROINTESTINAL DISORDERS",
            "DIZZINESS": "NERVOUS SYSTEM DISORDERS",
            "APPLICATION SITE PRURITUS": "GENERAL DISORDERS AND ADMINISTRATION SITE CONDITIONS",
            "RASH":     "SKIN AND SUBCUTANEOUS TISSUE DISORDERS",
            "FATIGUE":  "GENERAL DISORDERS AND ADMINISTRATION SITE CONDITIONS",
            "VOMITING": "GASTROINTESTINAL DISORDERS",
            "ERYTHEMA": "SKIN AND SUBCUTANEOUS TISSUE DISORDERS",
            "DIARRHOEA": "GASTROINTESTINAL DISORDERS",
            "PRURITUS": "SKIN AND SUBCUTANEOUS TISSUE DISORDERS",
        }
        return pd.DataFrame({
            "STUDYID":  "CDISCPILOT01",
            "USUBJID":  subjects,
            "ACTARM":   arms,
            "AETERM":   ae_terms,
            "AEDECOD":  ae_terms,
            "AESOC":    [soc_map[t] for t in ae_terms],
            "AESEV":    severities,
            "AESER":    rng.choice(["Y", "N"], n, p=[0.1, 0.9]),
            "TRTEMFL":  rng.choice(["Y", "N"], n, p=[0.85, 0.15]),
        })

    def _setup_openai_chain(self, api_key: str, model: str) -> None:
        """Configure a LangChain chain for OpenAI."""
        llm = ChatOpenAI(
            api_key     = api_key,
            model       = model,
            temperature = 0,          # Deterministic output for structured parsing
            max_tokens  = 200,
        )
        prompt = ChatPromptTemplate.from_messages([
            ("system", SYSTEM_PROMPT),
            ("human",  USER_PROMPT),
        ])
        self._chain = prompt | llm | StrOutputParser()

    def _call_llm(self, question: str) -> str:
        """Call the configured LLM and return raw string response."""
        if self._llm_mode == "openai" and self._chain is not None:
            return self._chain.invoke({
                "schema":   SCHEMA_STR,
                "question": question,
            })

        if self._llm_mode == "anthropic":
            client = anthropic_sdk.Anthropic(api_key=self._anthropic_key)
            msg = client.messages.create(
                model      = "claude-sonnet-4-20250514",
                max_tokens = 200,
                system     = SYSTEM_PROMPT.format(schema=SCHEMA_STR),
                messages   = [{"role": "user", "content": USER_PROMPT.format(question=question)}],
            )
            return msg.content[0].text

        # Mock mode
        result = _mock_llm_parse(question)
        return json.dumps(result)

    # ──────────────────────────────────────────────────────────────────────
    # Public API
    # ──────────────────────────────────────────────────────────────────────

    def parse_question(self, question: str) -> dict[str, str]:
        """
        Send a natural language question to the LLM and parse the JSON response.

        Returns
        -------
        dict with keys "target_column" and "filter_value"
        """
        raw_response = self._call_llm(question)

        # Extract JSON from response (handles cases where LLM adds markdown)
        json_match = re.search(r"\{.*?\}", raw_response, re.DOTALL)
        if not json_match:
            raise ValueError(f"LLM did not return valid JSON. Response was:\n{raw_response}")

        parsed = json.loads(json_match.group())

        # Validate required keys
        for key in ("target_column", "filter_value"):
            if key not in parsed:
                raise KeyError(f"LLM response missing required key: '{key}'")

        # Validate column exists in schema
        if parsed["target_column"] not in ADAE_SCHEMA:
            raise ValueError(
                f"LLM returned unknown column '{parsed['target_column']}'. "
                f"Valid columns: {list(ADAE_SCHEMA.keys())}"
            )

        return parsed

    def execute_filter(self, parsed_filter: dict[str, str]) -> pd.DataFrame:
        """
        Apply the parsed filter to the ADAE dataframe.

        Parameters
        ----------
        parsed_filter : dict with "target_column" and "filter_value"

        Returns
        -------
        Filtered DataFrame
        """
        col = parsed_filter["target_column"]
        val = parsed_filter["filter_value"]

        if col not in self.df.columns:
            raise ValueError(f"Column '{col}' not found in dataset. "
                             f"Available columns: {list(self.df.columns)}")

        # Case-insensitive string matching for robustness
        col_series = self.df[col].astype(str).str.upper()
        val_upper  = str(val).upper()

        return self.df[col_series == val_upper].copy()

    def query(self, question: str) -> QueryResult:
        """
        End-to-end: parse question -> execute filter -> return result.

        This is the main entry point for the agent.

        Parameters
        ----------
        question : str
            Natural language question about the AE dataset.

        Returns
        -------
        QueryResult dataclass with subject count, IDs, and filter details.
        """
        print(f"\n[QUERY] {question}")

        # Step 1: Prompt → LLM → JSON
        parsed_filter = self.parse_question(question)
        print(f"[PARSED] {parsed_filter}")

        # Step 2: JSON → Pandas filter
        filtered_df = self.execute_filter(parsed_filter)

        # Step 3: Aggregate results
        unique_subjects = filtered_df["USUBJID"].dropna().unique().tolist()
        n_subjects      = len(unique_subjects)

        return QueryResult(
            question      = question,
            parsed_filter = parsed_filter,
            n_subjects    = n_subjects,
            subject_ids   = sorted(unique_subjects),
            n_records     = len(filtered_df),
            llm_mode      = self._llm_mode,
        )


# =============================================================================
# 5. Test Script – 3 Example Queries
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print(" Clinical Trial Data Agent – Test Script")
    print("=" * 60)

    # Initialize agent
    # Uses API key from environment variable OPENAI_API_KEY if set,
    # otherwise falls back to mock mode
    agent = ClinicalTrialDataAgent(
        api_key           = os.getenv("OPENAI_API_KEY"),
        anthropic_api_key = os.getenv("ANTHROPIC_API_KEY"),
    )

    # ── Query 1: Severity filter ──────────────────────────────────────────
    result1 = agent.query(
        "Give me the subjects who had adverse events of Moderate severity."
    )
    print(result1)

    # ── Query 2: Body system filter ───────────────────────────────────────
    result2 = agent.query(
        "Which patients experienced adverse events related to the nervous system?"
    )
    print(result2)

    # ── Query 3: Specific AE term ─────────────────────────────────────────
    result3 = agent.query(
        "Find all subjects who reported nausea as an adverse event."
    )
    print(result3)

    print("\n=== All 3 queries completed successfully ===")
