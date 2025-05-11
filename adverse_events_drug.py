from typing import List
import requests
import pandas as pd
import time

from datasets import load_dataset
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# if getting bad gateway errors, try:
# Retry session setup
session = requests.Session()
retries = Retry(
    total=5,  # Retry up to 5 times
    backoff_factor=1,  # Wait 1s, then 2s, then 4s...
    status_forcelist=[502, 503, 504],
    allowed_methods=["GET"],
)
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_faers_adverse_events(drug_name: str, limit: int = 100) -> List:
    """
    Fetch and parse adverse event reports for a drug from openFDA FAERS.
    Adds serious outcome flags and event context fields.
    Args:
        drug_name (str): The name of the drug to search for.
        limit (int): The maximum number of records to fetch.
    Returns:
        List: A list of dictionaries containing adverse event data.
    """
    print(f"Fetching FAERS data for: {drug_name}")
    url = "https://api.fda.gov/drug/event.json"
    params = {"search": f'patient.drug.medicinalproduct:"{drug_name}"', "limit": limit}

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        results = response.json().get("results", [])

        # Extract detailed adverse event fields
        events = []
        for entry in results:
            patient = entry.get("patient", {})
            reactions = [r.get("reactionmeddrapt") for r in patient.get("reaction", [])]

            event = {
                "drug_name": drug_name,
                "safetyreportid": entry.get("safetyreportid"),
                "serious": entry.get("serious"),
                "seriousnessdeath": entry.get("seriousnessdeath"),
                "seriousnesshospitalization": entry.get("seriousnesshospitalization"),
                "seriousnessdisabling": entry.get("seriousnessdisabling"),
                "seriousnesslifethreatening": entry.get("seriousnesslifethreatening"),
                "seriousnesscongenitalanomali": entry.get("seriousnesscongenitalanomali"),
                "reactions": reactions,
                "receivedate": entry.get("receivedate"),
                "occurcountry": entry.get("occurcountry"),
            }
            events.append(event)
        print(f"Fetched {len(events)} events for {drug_name}")

        return events

    except requests.exceptions.RequestException as e:
        print(f"Failed to fetch data for {drug_name}: {e}")
        return []


def pull_faers_data(drug_list: List[str], output_filename_prefix: str) -> pd.DataFrame:
    """
    Pull FAERS data for a list of drugs and save to CSV.
    Args:
        drug_list (List[str]): List of drug names to fetch data for.
    """
    all_events = []
    for drug in drug_list:
        events = get_faers_adverse_events(drug, limit=100)
        all_events.extend(events)
        time.sleep(1)  # be polite to openFDA rate limits

    # Convert to DataFrame and show or save
    df = pd.DataFrame(all_events)
    df.to_csv(f"{output_filename_prefix}_faers_adverse_events.csv", index=False)
    print(df.head(), df.shape)
    return df


# example drug list
drug_list = ["pembrolizumab", "nivolumab", "ipilimumab"]

# test out fx
pull_faers_data(drug_list=drug_list, output_filename_prefix="test")


# load tahoe drug data metadata, gather the AEs for each unique drug
ds = load_dataset("tahoebio/Tahoe-100M", "drug_metadata")
drug_metadata = ds["train"].to_pandas()

# get adverse events for each unique drug in the tahoe dataset
tahoe_drug_list = drug_metadata["drug"].unique().tolist()
adverse_events = pull_faers_data(drug_list=tahoe_drug_list, output_filename_prefix="tahoe_drugs")
# adverse_events = pd.read_csv("tahoe_drugs_faers_adverse_events.csv")

# summarize the adverse events for each drug so we can try to predict the severity of a drug
def compute_severity_score(df: pd.DataFrame) -> float:
    weights = {
        "seriousnessdeath": 5,
        "seriousnesslifethreatening": 4,
        "seriousnesshospitalization": 3,
        "seriousnessdisabling": 2,
        "seriousnesscongenitalanomali": 1,
    }
    score = 0
    for k, w in weights.items():
        score += df[k].fillna("0").astype(int).sum() * w
    return score


def compute_seriousness_ratio(df: pd.DataFrame) -> float:
    """
    Compute the ratio of serious to non-serious adverse events.
    """
    seriousness_ae_ratio = (
        df.groupby("drug_name")["serious"]
        .apply(lambda x: x.fillna("0").astype(int).mean())
        .reset_index(name="mean_serious_flag")
    )
    return seriousness_ae_ratio


ae_score_df = adverse_events.groupby("drug_name").apply(compute_severity_score).reset_index(name="ae_severity_score")
seriousness_ratio = compute_seriousness_ratio(adverse_events)

# write out outputs
ae_score_df.to_csv("tahoe_drug_ae_severity_score.csv", index=False)
seriousness_ratio.to_csv("tahoe_drug_ae_seriousness_ratio.csv", index=False)
