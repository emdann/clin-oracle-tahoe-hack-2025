import os
import requests
import pandas as pd
import json
import time
import re
from typing import List, Dict, Any, Union, Optional

import anthropic

class DrugTrialExtractor:
    """
    A class to extract information about which drugs are in clinical trials and their indications
    using PubChem IDs and clinical trials data sources.
    """
    
    def __init__(self):
        # Base URLs for APIs
        self.pubchem_base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.pubchem_view_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound"
        self.clinicaltrials_api_url = "https://clinicaltrials.gov/api/v2/studies"
        
        # For rate limiting
        self.request_delay = 0.5  # seconds between API requests
    
    def get_drug_info_from_pubchem(self, pubchem_id: str) -> Dict[str, Any]:
        """
        Get basic drug information from PubChem using a PubChem ID (CID)
        
        Parameters:
        -----------
        pubchem_id : str
            The PubChem Compound ID (CID)
            
        Returns:
        --------
        Dict[str, Any]
            Dictionary with drug information
        """
        try:
            # Validate the PubChem ID format
            try:
                # Ensure it's numeric and can be converted to integer
                int(pubchem_id)
            except ValueError:
                return {
                    "pubchem_id": pubchem_id, 
                    "error": "Invalid PubChem ID format - must be numeric"
                }
            
            # First check if the compound exists with a simple request
            check_url = f"{self.pubchem_base_url}/compound/cid/{pubchem_id}/cids/JSON"
            check_response = requests.get(check_url)
            
            if check_response.status_code != 200:
                error_msg = f"Error code {check_response.status_code}"
                if check_response.status_code == 400:
                    error_msg = "PubChem ID doesn't exist or has been deprecated"
                elif check_response.status_code == 404:
                    error_msg = "PubChem ID not found"
                elif check_response.status_code == 429:
                    error_msg = "Rate limit exceeded - try again later"
                    
                return {"pubchem_id": pubchem_id, "error": error_msg}
            
            # Get drug properties (name, formula, etc.)
            property_url = f"{self.pubchem_base_url}/compound/cid/{pubchem_id}/property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES/JSON"
            property_response = requests.get(property_url)
            
            if property_response.status_code != 200:
                print(f"Error retrieving properties for PubChem ID {pubchem_id}: {property_response.status_code}")
                print(f"Response: {property_response.text}")
                return {"pubchem_id": pubchem_id, "error": "Failed to retrieve properties"}
            
            property_data = property_response.json()
            
            # Add delay to respect rate limits
            time.sleep(self.request_delay)
            
            # Get synonyms to find common drug names
            synonyms_url = f"{self.pubchem_base_url}/compound/cid/{pubchem_id}/synonyms/JSON"
            synonyms_response = requests.get(synonyms_url)
            
            if synonyms_response.status_code != 200:
                print(f"Error retrieving synonyms for PubChem ID {pubchem_id}: {synonyms_response.status_code}")
                synonyms = []
            else:
                synonyms_data = synonyms_response.json()
                synonyms = synonyms_data.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
            
            # Add delay to respect rate limits
            time.sleep(self.request_delay)
            
            # Compile the drug information
            props = property_data.get("PropertyTable", {}).get("Properties", [{}])[0]
            
            # Find the common name (usually a drug name) from synonyms
            # This is heuristic-based and prioritizes known drug name patterns
            common_name = None
            drug_name_candidates = []
            
            if synonyms:
                for syn in synonyms:
                    # Skip very long names and suspected CAS numbers or identifiers
                    if len(syn) > 50 or bool(re.match(r'^\d+-\d+-\d+$', syn)):
                        continue
                        
                    # Prioritize names that look like drug names (no special chars, not too long)
                    if len(syn) < 20 and syn.isalpha() and syn[0].isupper():
                        drug_name_candidates.append(syn)
                
                # Sort by length (shorter names first)
                drug_name_candidates.sort(key=len)
                
                if drug_name_candidates:
                    common_name = drug_name_candidates[0]
                else:
                    # Fallback to any reasonable synonym
                    sorted_synonyms = sorted([s for s in synonyms if len(s) < 30], key=len)
                    if sorted_synonyms:
                        common_name = sorted_synonyms[0]
            
            drug_info = {
                "pubchem_id": pubchem_id,
                "iupac_name": props.get("IUPACName", ""),
                "common_name": common_name,
                "molecular_formula": props.get("MolecularFormula", ""),
                "molecular_weight": props.get("MolecularWeight", ""),
                "canonical_smiles": props.get("CanonicalSMILES", ""),
                "synonyms": synonyms[:10] if len(synonyms) > 10 else synonyms  # Limit to 10 synonyms
            }
            
            return drug_info
            
        except Exception as e:
            print(f"Error in get_drug_info_from_pubchem for ID {pubchem_id}: {str(e)}")
            return {"pubchem_id": pubchem_id, "error": str(e)}
    
    def get_clinical_trials_sections(self, pubchem_id: str) -> Optional[Dict[str, Any]]:
        """
        Get the clinical trials sections from PubChem for a given compound
        
        Parameters:
        -----------
        pubchem_id : str
            The PubChem Compound ID (CID)
            
        Returns:
        --------
        Optional[Dict[str, Any]]
            Dictionary with clinical trials sections or None if not found
        """
        try:
            # Using PUG-View to get the clinical trials sections
            url = f"{self.pubchem_view_url}/{pubchem_id}/JSON?response_type=display"
            response = requests.get(url)
            
            if response.status_code != 200:
                print(f"Error retrieving clinical trials sections for PubChem ID {pubchem_id}: {response.status_code}")
                return None
            
            data = response.json()
            
            # Navigate through the JSON to find the clinical trials section
            record = data.get("Record", {})
            section = record.get("Section", [])
            
            clinical_trials_data = None
            
            # Look for "Drug and Medication Information" section
            for s in section:
                if s.get("TOCHeading") == "Drug and Medication Information":
                    drug_section = s.get("Section", [])
                    for ds in drug_section:
                        if ds.get("TOCHeading") == "Clinical Trials":
                            clinical_trials_data = ds
                            break
                    break
            
            return clinical_trials_data
            
        except Exception as e:
            print(f"Error in get_clinical_trials_sections for ID {pubchem_id}: {str(e)}")
            return None
    
    def extract_trials_from_pubchem(self, pubchem_id: str) -> List[Dict[str, Any]]:
        """
        Extract clinical trials information from PubChem for a given compound ID
        
        Parameters:
        -----------
        pubchem_id : str
            The PubChem Compound ID (CID)
            
        Returns:
        --------
        List[Dict[str, Any]]
            List of dictionaries with clinical trial information
        """
        clinical_trials_section = self.get_clinical_trials_sections(pubchem_id)
        
        if not clinical_trials_section:
            return []
        
        trial_list = []
        
        # Process the clinical trials section to extract trials information
        try:
            section_list = clinical_trials_section.get("Section", [])
            
            for section in section_list:
                source_name = section.get("TOCHeading", "Unknown Source")
                
                # Look for information about trials in this source
                info_list = section.get("Information", [])
                
                for info in info_list:
                    if "ExternalTableName" in info.get("Value", {}):
                        # This means there are trials in this source
                        trial_list.append({
                            "pubchem_id": pubchem_id,
                            "source": source_name,
                            "has_trials": True,
                            "trial_count": info.get("Value", {}).get("ExternalTableNumRows", "Unknown")
                        })
        
        except Exception as e:
            print(f"Error extracting trials from PubChem for ID {pubchem_id}: {str(e)}")
        
        return trial_list
        
    def search_clinicaltrials_gov(self, drug_name: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """
        Search ClinicalTrials.gov for trials involving a specific drug
        
        Parameters:
        -----------
        drug_name : str
            The name of the drug to search for
        max_results : int, optional
            Maximum number of results to return (default 100)
            
        Returns:
        --------
        List[Dict[str, Any]]
            List of dictionaries with clinical trial information
        """
        try:
            # Skip names that are too long or complex as they'll likely cause a 400 error
            if len(drug_name) > 150 or drug_name.count('[') > 2:
                print(f"Skipping search for complex name: {drug_name}")
                return []
                
            # Using the new ClinicalTrials.gov API v2
            params = {
                "query.term": drug_name,
                "pageSize": min(max_results, 100)  # API limits to 1000 per request
            }
            
            trials = []
            
            # Make initial request
            response = requests.get(self.clinicaltrials_api_url, params=params)
            
            if response.status_code != 200:
                print(f"Error searching ClinicalTrials.gov for {drug_name}: {response.status_code}")
                if response.status_code == 400:
                    # If it's a 400 Bad Request, the drug name is likely invalid for the API
                    print(f"Drug name '{drug_name}' is not valid for the ClinicalTrials.gov API.")
                return []
            
            data = response.json()
            studies = data.get('studies', [])
            
            # Process all studies from the first page
            for study in studies:
                try:
                    # Extract the protocol section which contains most of the important information
                    protocol = study.get("protocolSection", {})
                    
                    # Extract identification information
                    identification = protocol.get("identificationModule", {})
                    nct_id = identification.get("nctId", "Unknown")
                    brief_title = identification.get("briefTitle", "Unknown")
                    
                    # Extract status information
                    status_module = protocol.get("statusModule", {})
                    overall_status = status_module.get("overallStatus", "Unknown")
                    
                    # FIXED: Extract phase information from designModule.phases array
                    design_module = protocol.get("designModule", {})
                    phases = design_module.get("phases", [])
                    phase = phases[0] if phases else "Unknown"
                    
                    # Extract conditions
                    conditions_module = protocol.get("conditionsModule", {})
                    conditions = conditions_module.get("conditions", [])
                    
                    # Extract interventions (drugs, etc.)
                    interventions = []
                    intervention_module = protocol.get("armsInterventionsModule", {})
                    intervention_list = intervention_module.get("interventions", [])
                    
                    for intervention in intervention_list:
                        intervention_name = intervention.get("name", "")
                        intervention_type = intervention.get("type", "")
                        intervention_description = intervention.get("description", "")
                        
                        interventions.append({
                            "name": intervention_name,
                            "type": intervention_type,
                            "description": intervention_description
                        })
                    
                    # Create a dictionary for this trial
                    trial_info = {
                        "nct_id": nct_id,
                        "title": brief_title,
                        "status": overall_status,
                        "phase": phase,
                        "conditions": conditions,
                        "interventions": interventions
                    }
                    
                    trials.append(trial_info)
                    
                except Exception as e:
                    print(f"Error processing study {study.get('protocolSection', {}).get('identificationModule', {}).get('nctId', 'Unknown')}: {str(e)}")
            
            # Check if there are more pages
            next_page_token = data.get("nextPageToken")
            
            # Continue getting data if there's a next page and we haven't reached max_results
            while next_page_token and len(trials) < max_results:
                # Add delay to respect rate limits
                time.sleep(self.request_delay)
                
                # Update parameters with the next page token
                params["pageToken"] = next_page_token
                
                # Make the request for the next page
                response = requests.get(self.clinicaltrials_api_url, params=params)
                
                if response.status_code != 200:
                    print(f"Error retrieving next page for {drug_name}: {response.status_code}")
                    break
                
                data = response.json()
                studies = data.get("studies", [])
                
                # Process all studies from this page
                for study in studies:
                    if len(trials) >= max_results:
                        break
                        
                    try:
                        # Extract the protocol section
                        protocol = study.get("protocolSection", {})
                        
                        # Extract identification information
                        identification = protocol.get("identificationModule", {})
                        nct_id = identification.get("nctId", "Unknown")
                        brief_title = identification.get("briefTitle", "Unknown")
                        
                        # Extract status information
                        status_module = protocol.get("statusModule", {})
                        overall_status = status_module.get("overallStatus", "Unknown")
                        
                        # FIXED: Extract phase information from designModule.phases array
                        design_module = protocol.get("designModule", {})
                        phases = design_module.get("phases", [])
                        phase = phases[0] if phases else "Unknown"
                        
                        # Extract conditions
                        conditions_module = protocol.get("conditionsModule", {})
                        conditions = conditions_module.get("conditions", [])
                        
                        # Extract interventions (drugs, etc.)
                        interventions = []
                        intervention_module = protocol.get("armsInterventionsModule", {})
                        intervention_list = intervention_module.get("interventions", [])
                        
                        for intervention in intervention_list:
                            intervention_name = intervention.get("name", "")
                            intervention_type = intervention.get("type", "")
                            intervention_description = intervention.get("description", "")
                            
                            interventions.append({
                                "name": intervention_name,
                                "type": intervention_type,
                                "description": intervention_description
                            })
                        
                        # Create a dictionary for this trial
                        trial_info = {
                            "nct_id": nct_id,
                            "title": brief_title,
                            "status": overall_status,
                            "phase": phase,
                            "conditions": conditions,
                            "interventions": interventions
                        }
                        
                        trials.append(trial_info)
                        
                    except Exception as e:
                        print(f"Error processing study {study.get('protocolSection', {}).get('identificationModule', {}).get('nctId', 'Unknown')}: {str(e)}")
                
                # Update the next page token for the next iteration
                next_page_token = data.get("nextPageToken")
            
            return trials
            
        except Exception as e:
            print(f"Error in search_clinicaltrials_gov for drug {drug_name}: {str(e)}")
            return []
    
    def get_trials_for_drug(self, pubchem_id: str, max_results: int = 100) -> Dict[str, Any]:
        """
        Get complete clinical trial information for a drug with the given PubChem ID
        
        Parameters:
        -----------
        pubchem_id : str
            The PubChem Compound ID (CID)
        max_results : int, optional
            Maximum number of clinical trials to retrieve per drug (default 100)
            
        Returns:
        --------
        Dict[str, Any]
            Dictionary with drug information and associated clinical trials
        """
        # Get drug information
        drug_info = self.get_drug_info_from_pubchem(pubchem_id)
        
        # Check if PubChem has clinical trials information
        pubchem_trials_info = self.extract_trials_from_pubchem(pubchem_id)
        
        # Search ClinicalTrials.gov using drug names
        trials = []
        
        # Try with common name first as it's most likely to work
        if drug_info.get("common_name"):
            print(f"Searching ClinicalTrials.gov for common name: {drug_info['common_name']}")
            trials.extend(self.search_clinicaltrials_gov(drug_info["common_name"], max_results))
        
        # Try with a filtered list of synonyms (prioritizing drug names over chemical identifiers)
        if len(trials) < max_results and drug_info.get("synonyms"):
            drug_synonyms = []
            for synonym in drug_info.get("synonyms", []):
                # Skip long and complex names, molecular identifiers, and registry numbers
                if (len(synonym) < 30 and 
                    not any(char in synonym for char in "[](){}-=") and
                    not synonym.isdigit() and
                    not bool(re.match(r'^\d+-\d+-\d+$', synonym))):  # Skip CAS registry numbers
                    drug_synonyms.append(synonym)
            
            # Prioritize shorter names as they're more likely to be common drug names
            drug_synonyms.sort(key=len)
            
            # Try up to 3 prioritized synonyms
            for synonym in drug_synonyms[:3]:
                if len(trials) >= max_results:
                    break
                
                print(f"Searching ClinicalTrials.gov for synonym: {synonym}")
                try:
                    more_trials = self.search_clinicaltrials_gov(synonym, max_results - len(trials))
                    
                    # Filter to avoid duplicates
                    existing_nct_ids = {t["nct_id"] for t in trials}
                    for trial in more_trials:
                        if trial["nct_id"] not in existing_nct_ids:
                            trials.append(trial)
                            existing_nct_ids.add(trial["nct_id"])
                except Exception as e:
                    print(f"Error searching for synonym {synonym}: {str(e)}")
        
        # As a last resort, try with IUPAC name (least likely to work with clinical trials API)
        if len(trials) < max_results and len(trials) == 0 and drug_info.get("iupac_name"):
            # Only use IUPAC name if it's reasonably short
            if len(drug_info["iupac_name"]) < 100:  # Avoid extremely long IUPAC names
                print(f"Searching ClinicalTrials.gov for IUPAC name: {drug_info['iupac_name']}")
                try:
                    more_trials = self.search_clinicaltrials_gov(drug_info["iupac_name"], max_results - len(trials))
                    
                    # Filter to avoid duplicates
                    existing_nct_ids = {t["nct_id"] for t in trials}
                    for trial in more_trials:
                        if trial["nct_id"] not in existing_nct_ids:
                            trials.append(trial)
                            existing_nct_ids.add(trial["nct_id"])
                except Exception as e:
                    print(f"Error searching for IUPAC name: {str(e)}")
        
        # Compile results
        result = {
            "drug_info": drug_info,
            "pubchem_trials_info": pubchem_trials_info,
            "clinicaltrials_gov_trials": trials,
        }
        
        return result
    
    def process_multiple_drugs(self, pubchem_ids: List[str], max_results_per_drug: int = 100) -> Dict[str, List[Dict[str, Any]]]:
        """
        Process multiple drugs and get their clinical trial information
        
        Parameters:
        -----------
        pubchem_ids : List[str]
            List of PubChem Compound IDs (CIDs)
        max_results_per_drug : int, optional
            Maximum number of clinical trials to retrieve per drug (default 100)
            
        Returns:
        --------
        Dict[str, List[Dict[str, Any]]]
            Dictionary with results for each drug
        """
        results = {}
        
        for pubchem_id in pubchem_ids:
            print(f"Processing PubChem ID: {pubchem_id}")
            drug_results = self.get_trials_for_drug(pubchem_id, max_results_per_drug)
            results[pubchem_id] = drug_results
            
            # Add delay to respect rate limits
            time.sleep(self.request_delay)
        
        return results
    
    def save_results_to_csv(self, results: Dict[str, Any], output_prefix: str = "drug_trials") -> Dict[str, str]:
        """
        Save the results to CSV files
        
        Parameters:
        -----------
        results : Dict[str, Any]
            Results from process_multiple_drugs
        output_prefix : str, optional
            Prefix for output CSV files (default "drug_trials")
            
        Returns:
        --------
        Dict[str, str]
            Dictionary with paths to the output files
        """
        # Create DataFrames
        drug_info_rows = []
        all_trials = []
        pubchem_trials_rows = []
        
        for pubchem_id, drug_data in results.items():
            # Add drug info
            drug_info = drug_data.get("drug_info", {})
            drug_info_rows.append(drug_info)
            
            # Add PubChem trials info
            for trial_info in drug_data.get("pubchem_trials_info", []):
                pubchem_trials_rows.append(trial_info)
            
            # Add ClinicalTrials.gov trials
            for trial in drug_data.get("clinicaltrials_gov_trials", []):
                trial_copy = trial.copy()
                
                # Add drug info to the trial
                trial_copy["pubchem_id"] = pubchem_id
                trial_copy["drug_name"] = drug_info.get("common_name", drug_info.get("iupac_name", ""))
                
                # Convert lists to strings for CSV
                if "conditions" in trial_copy:
                    trial_copy["conditions"] = "; ".join(trial_copy["conditions"])
                
                if "interventions" in trial_copy:
                    interventions_list = trial_copy["interventions"]
                    intervention_strings = []
                    
                    for intervention in interventions_list:
                        int_str = f"{intervention.get('name', '')} ({intervention.get('type', '')})"
                        intervention_strings.append(int_str)
                    
                    trial_copy["interventions"] = "; ".join(intervention_strings)
                
                all_trials.append(trial_copy)
        
        # Create DataFrames
        drug_info_df = pd.DataFrame(drug_info_rows)
        trials_df = pd.DataFrame(all_trials)
        pubchem_trials_df = pd.DataFrame(pubchem_trials_rows)
        
        # Save to CSV
        drug_info_path = f"{output_prefix}_drug_info.csv"
        trials_path = f"{output_prefix}_clinical_trials.csv"
        pubchem_trials_path = f"{output_prefix}_pubchem_trials_info.csv"
        
        drug_info_df.to_csv(drug_info_path, index=False)
        trials_df.to_csv(trials_path, index=False)
        pubchem_trials_df.to_csv(pubchem_trials_path, index=False)
        
        return {
            "drug_info": drug_info_path,
            "clinical_trials": trials_path,
            "pubchem_trials_info": pubchem_trials_path
        }

# In your shared code:
def get_anthropic_client():
    # Get API key from environment variable
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    
    if not api_key:
        raise ValueError("Missing ANTHROPIC_API_KEY environment variable")
        
    return anthropic.Anthropic(api_key=api_key)

def standardize_medical_conditions(disease_list):
    """
    Standardizes a list of medical conditions using Claude API to identify and group similar conditions.
    
    Args:
        disease_list (list): List of medical condition strings to standardize
        api_key (str): Your Anthropic API key
        
    Returns:
        pandas.DataFrame: DataFrame with original conditions and their standardized labels
    """
    # Initialize the Claude client
    client = get_anthropic_client()
    
    # Create chunks if your list is very large (Claude has context limits)
    chunk_size = 50  # Adjust based on your needs
    all_groups = {}
    
    # Remove duplicates to reduce API costs while processing
    unique_diseases = list(set(disease_list))
    
    for i in range(0, len(unique_diseases), chunk_size):
        chunk = unique_diseases[i:i+chunk_size]
        
        # Format the prompt for Claude
        formatted_diseases = "\n".join([f"- {d}" for d in chunk])
        
        prompt = f"""
        Here's a list of medical conditions:
        {formatted_diseases}
        
        Please group these conditions into standardized categories where entries refer to the same basic condition. 
        For each group, select the most appropriate, specific, and concise label.
        
        IMPORTANT: Only group conditions when there's a clear case for doing so. When a condition is unique or 
        doesn't clearly fit with others, keep it as its own separate category.
        
        Examples:
        - "Prostate Cancer" and "Prostate Cancer; Prostate Adenocarcinoma" can be grouped as "Prostate Cancer"
        - But "Small Cell Lung Cancer" should NOT be grouped with "Non-Small Cell Lung Cancer" as they are distinct conditions
        - "Diabetes" and "Diabetes Mellitus Type 2" can be grouped, but should use the more specific "Diabetes Mellitus Type 2" as the label
        
        Format your response as a JSON dictionary where keys are the standardized labels and values are 
        lists of all original terms that should map to that label. Include every term from the input list.
        """
        
        # Call Claude API
        message = client.messages.create(
            model="claude-3-7-sonnet-20250219",
            max_tokens=4000,
            temperature=0,  # Keep it deterministic
            system="You are a medical terminology expert. Follow instructions exactly.",
            messages=[
                {"role": "user", "content": prompt}
            ]
        )
        
        # Extract JSON from response
        response_content = message.content[0].text
        json_match = re.search(r'```(?:json)?\s*([\s\S]*?)\s*```', response_content)
        
        if json_match:
            json_str = json_match.group(1)
        else:
            # If no code block, try to find JSON directly
            json_str = response_content
        
        try:
            chunk_groups = json.loads(json_str)
            all_groups.update(chunk_groups)
        except json.JSONDecodeError:
            print(f"Warning: Could not parse JSON for chunk {i}. Skipping this chunk.")
            continue
    
    # Convert to a mapping dictionary
    condition_mapping = {}
    for standard_label, variants in all_groups.items():
        for variant in variants:
            condition_mapping[variant] = standard_label
    
    # Apply mapping to original list (preserving order and duplicates)
    mapped_diseases = [condition_mapping.get(disease, disease) for disease in disease_list]
    
    # Create a DataFrame with the results
    result_df = pd.DataFrame({
        "Original": disease_list,
        "Standardized": mapped_diseases
    })
    
    return result_df
