import pandas as pd
import re

def load_indications(file_path):
    """
    Load indications from a TSV file.
    """
    try:
        # Read the TSV file
        df = pd.read_csv(file_path, sep='\t', na_values=['', 'NA', 'N/A'])
        return df
    except Exception as e:
        print(f"Error loading file: {e}")
        return None

def map_indication_to_organ(indication, areas=None):
    """
    Map a medical indication to an organ.
    Returns the organ name if a confident mapping exists, otherwise None.
    """
    # List of target organs
    target_organs = [
        'Bowel', 'Lung', 'Esophagus/Stomach', 'Pancreas', 'Skin', 
        'Uterus', 'Breast', 'Ovary/Fallopian Tube', 'Cervix', 
        'CNS/Brain', 'Liver', 'Kidney', 'Peripheral Nervous System', 
        'Vulva/Vagina', 'Bladder/Urinary Tract'
    ]
    
    # Dictionary mapping indication keywords to organs
    keyword_to_organ = {
        # Bowel related
        'bowel': 'Bowel',
        'intestinal': 'Bowel',
        'colon': 'Bowel',
        'colorectal': 'Bowel',
        'rectal': 'Bowel',
        'anus': 'Bowel',
        'colitis': 'Bowel',
        'crohn': 'Bowel',
        'inflammatory bowel': 'Bowel',
        'diverticulitis': 'Bowel',
        'ileus': 'Bowel',
        'short bowel': 'Bowel',
        
        # Lung related
        'lung': 'Lung',
        'pulmonary': 'Lung',
        'respiratory': 'Lung',
        'bronch': 'Lung',
        'pneumonia': 'Lung',
        'emphysema': 'Lung',
        'asthma': 'Lung',
        'chronic obstructive pulmonary': 'Lung',
        'idiopathic pulmonary fibrosis': 'Lung',
        
        # Esophagus/Stomach related
        'esophag': 'Esophagus/Stomach',
        'stomach': 'Esophagus/Stomach',
        'gastric': 'Esophagus/Stomach',
        'gastroesophageal': 'Esophagus/Stomach',
        'duodenal': 'Esophagus/Stomach',
        'peptic ulcer': 'Esophagus/Stomach',
        'gastritis': 'Esophagus/Stomach',
        'barrett': 'Esophagus/Stomach',
        
        # Pancreas related
        'pancrea': 'Pancreas',
        'exocrine pancreatic': 'Pancreas',
        
        # Skin related
        'skin': 'Skin',
        'dermat': 'Skin',
        'cutaneous': 'Skin',
        'melanoma': 'Skin',
        'psoriasis': 'Skin',
        'eczema': 'Skin',
        'acne': 'Skin',
        'rosacea': 'Skin',
        'alopecia': 'Skin',
        'vitiligo': 'Skin',
        'keratosis': 'Skin',
        'keloid': 'Skin',
        'ichthyosis': 'Skin',
        'pruritus': 'Skin',
        'urticaria': 'Skin',
        'pemphigus': 'Skin',
        'pemphigoid': 'Skin',
        'scleroderma': 'Skin',
        
        # Uterus related
        'uter': 'Uterus',
        'endometri': 'Uterus',
        'myometri': 'Uterus',
        'leiomyoma': 'Uterus',
        
        # Breast related
        'breast': 'Breast',
        'mammary': 'Breast',
        
        # Ovary/Fallopian tube related
        'ovar': 'Ovary/Fallopian Tube',
        'fallopian': 'Ovary/Fallopian Tube',
        'polycystic ovary': 'Ovary/Fallopian Tube',
        
        # Cervix related
        'cervix': 'Cervix',
        'cervical': 'Cervix',
        
        # CNS/Brain related
        'brain': 'CNS/Brain',
        'central nervous': 'CNS/Brain',
        'cns': 'CNS/Brain',
        'cerebr': 'CNS/Brain',
        'alzheimer': 'CNS/Brain',
        'parkinson': 'CNS/Brain',
        'glioma': 'CNS/Brain',
        'glioblastoma': 'CNS/Brain',
        'medulloblastoma': 'CNS/Brain',
        'meningioma': 'CNS/Brain',
        'epilepsy': 'CNS/Brain',
        'dementia': 'CNS/Brain',
        'encephalitis': 'CNS/Brain',
        'encephalopathy': 'CNS/Brain',
        'huntington': 'CNS/Brain',
        'multiple sclerosis': 'CNS/Brain',
        'stroke': 'CNS/Brain',
        
        # Liver related
        'liver': 'Liver',
        'hepat': 'Liver',
        'biliary': 'Liver',
        'cholang': 'Liver',
        'cirrhosis': 'Liver',
        'fatty liver': 'Liver',
        
        # Kidney related
        'kidney': 'Kidney',
        'renal': 'Kidney',
        'nephro': 'Kidney',
        'nephri': 'Kidney',
        'glomerulo': 'Kidney',
        'polycystic kidney': 'Kidney',
        
        # Peripheral Nervous System related
        'peripheral nervous': 'Peripheral Nervous System',
        'neuropathy': 'Peripheral Nervous System',
        'neuralgia': 'Peripheral Nervous System',
        
        # Vulva/Vagina related
        'vulva': 'Vulva/Vagina',
        'vagina': 'Vulva/Vagina',
        'vulvar': 'Vulva/Vagina',
        'vaginal': 'Vulva/Vagina',
        'vaginitis': 'Vulva/Vagina',
        
        # Bladder/Urinary Tract related
        'bladder': 'Bladder/Urinary Tract',
        'urinary': 'Bladder/Urinary Tract',
        'urethral': 'Bladder/Urinary Tract',
        'urethra': 'Bladder/Urinary Tract',
        'cystitis': 'Bladder/Urinary Tract',
    }
    
    # Indications that are too general or affect multiple systems
    too_general = [
        'Neoplasms', 'Inflammation', 'Pain', 'Infection', 'Fever', 
        'Autoimmune Diseases', 'Immune System Diseases', 
        'Nervous System Diseases', 'Metabolic Diseases', 
        'Arthritis', 'Diabetes Mellitus', 'Hypertension',
        'Depression', 'Anxiety Disorders', 'HIV Infections',
        'Sepsis', 'Wound Healing', 'Cachexia', 'Fatigue',
        'Nausea', 'Vomiting', 'Multiple Myeloma', 'Leukemia', 'Lymphoma'
    ]
    
    # Specific disease to organ mapping
    specific_disease_mapping = {
        'Ulcerative Colitis': 'Bowel',
        'Crohn Disease': 'Bowel',
        'Asthma': 'Lung',
        'Psoriasis': 'Skin',
        'Lupus Erythematosus, Cutaneous': 'Skin',
        'Lupus Nephritis': 'Kidney',
        'Alzheimer Disease': 'CNS/Brain',
        'Parkinson Disease': 'CNS/Brain',
        'Amyotrophic Lateral Sclerosis': 'CNS/Brain',
        'Huntington Disease': 'CNS/Brain',
        'Non-alcoholic Fatty Liver Disease': 'Liver',
        'Carcinoma, Non-Small-Cell Lung': 'Lung',
        'Small Cell Lung Carcinoma': 'Lung',
        'Liver Cirrhosis': 'Liver',
        'Liver Neoplasms': 'Liver',
        'Liver Cirrhosis, Biliary': 'Liver',
        'Breast Neoplasms': 'Breast',
        'Ovarian Neoplasms': 'Ovary/Fallopian Tube',
        'Brain Neoplasms': 'CNS/Brain',
        'Brain Injuries': 'CNS/Brain',
        'Brain Injuries, Traumatic': 'CNS/Brain',
        'Brain Ischemia': 'CNS/Brain',
        'Pancreatic Neoplasms': 'Pancreas',
        'Urinary Bladder Neoplasms': 'Bladder/Urinary Tract',
        'Stomach Neoplasms': 'Esophagus/Stomach',
        'Uterine Cervical Neoplasms': 'Cervix',
        'Endometrial Neoplasms': 'Uterus',
        'Kidney Neoplasms': 'Kidney',
        'Colorectal Neoplasms': 'Bowel',
        'Inflammatory Bowel Diseases': 'Bowel',
        'COVID-19': 'Lung',
        'Severe Acute Respiratory Syndrome': 'Lung',
        'Pulmonary Disease, Chronic Obstructive': 'Lung',
        'Idiopathic Pulmonary Fibrosis': 'Lung',
        'Respiratory Distress Syndrome': 'Lung',
        'Esophageal Neoplasms': 'Esophagus/Stomach',
        'Acne Vulgaris': 'Skin',
        'Psoriatic Arthritis': 'Skin',
        'Alopecia': 'Skin',
        'Alopecia Areata': 'Skin',
        'Vitiligo': 'Skin',
        'Rosacea': 'Skin',
        'Hidradenitis Suppurativa': 'Skin',
        'Pemphigus': 'Skin',
        'Pemphigoid, Bullous': 'Skin',
        'Scleroderma, Systemic': 'Skin',
        'Scleroderma, Diffuse': 'Skin',
        'Scleroderma, Limited': 'Skin',
        'Fallopian Tube Neoplasms': 'Ovary/Fallopian Tube',
        'Polycystic Ovary Syndrome': 'Ovary/Fallopian Tube',
        'Endometriosis': 'Uterus',
        'Leiomyoma': 'Uterus',
        'Carcinoma, Squamous Cell': 'Skin',
        'Neoplasms, Basal Cell': 'Skin',
        'Melanoma': 'Skin',
        'Glioblastoma': 'CNS/Brain',
        'Medulloblastoma': 'CNS/Brain',
        'Irritable Bowel Syndrome': 'Bowel',
        'Diabetic Nephropathies': 'Kidney',
        'Diabetic Neuropathies': 'Peripheral Nervous System',
        'Multiple Sclerosis': 'CNS/Brain',
    }
    
    # First, check if the indication is too general
    if indication in too_general:
        return None
        
    # Check if the indication is in our specific disease mapping
    if indication in specific_disease_mapping:
        return specific_disease_mapping[indication]
    
    # Check for keyword matches
    indication_lower = indication.lower()
    for keyword, organ in keyword_to_organ.items():
        if keyword in indication_lower:
            return organ
            
    # Use the 'areas' context for additional mapping
    if areas:
        areas_lower = str(areas).lower()
        if 'digestive' in areas_lower:
            if any(term in indication_lower for term in ['bowel', 'colon', 'rectal', 'intestin']):
                return 'Bowel'
            elif any(term in indication_lower for term in ['liver', 'hepat', 'biliary']):
                return 'Liver'
            elif 'pancrea' in indication_lower:
                return 'Pancreas'
            elif any(term in indication_lower for term in ['stomach', 'gastric', 'esophag']):
                return 'Esophagus/Stomach'
        elif 'respiratory' in areas_lower:
            return 'Lung'
        elif 'neurology' in areas_lower:
            if 'peripheral' in indication_lower:
                return 'Peripheral Nervous System'
            else:
                return 'CNS/Brain'
        elif 'dermatology' in areas_lower:
            return 'Skin'
        elif 'oncology' in areas_lower:
            # For oncology, need more specific organ information
            if any(term in indication_lower for term in ['lung', 'pulmonary', 'bronch']):
                return 'Lung'
            elif any(term in indication_lower for term in ['brain', 'glioma', 'cerebr']):
                return 'CNS/Brain'
            elif any(term in indication_lower for term in ['liver', 'hepatic']):
                return 'Liver'
            elif any(term in indication_lower for term in ['skin', 'melanoma', 'basal cell']):
                return 'Skin'
            elif any(term in indication_lower for term in ['breast', 'mammary']):
                return 'Breast'
            elif any(term in indication_lower for term in ['colon', 'rectal', 'colorectal']):
                return 'Bowel'
            elif any(term in indication_lower for term in ['pancreas', 'pancreatic']):
                return 'Pancreas'
            elif any(term in indication_lower for term in ['kidney', 'renal']):
                return 'Kidney'
            elif any(term in indication_lower for term in ['bladder', 'urinary']):
                return 'Bladder/Urinary Tract'
            elif any(term in indication_lower for term in ['cervix', 'cervical']):
                return 'Cervix'
            elif any(term in indication_lower for term in ['uterus', 'uterine', 'endometri']):
                return 'Uterus'
            elif any(term in indication_lower for term in ['ovary', 'ovarian', 'fallopian']):
                return 'Ovary/Fallopian Tube'
            
    # If no confident mapping is found, return None
    return None

def main():
    # Load indications from TSV file
    file_path = '../genetic_support/data/indic.tsv'  # Update this path if necessary
    indications_df = load_indications(file_path)
    indications_df = indications_df.dropna(subset=['indication_mesh_term'])
    
    if indications_df is not None:
        # Add a new column for organ mapping
        indications_df['organ'] = indications_df.apply(
            lambda row: map_indication_to_organ(row['indication_mesh_term'], row.get('areas')), 
            axis=1
        )
        
        # Display indications with their mapped organs
        mapped_df = indications_df[['indication_mesh_term', 'organ']]
        all_mappings = mapped_df.copy()  # Save all mappings for reference
        
        # Show only confident mappings
        mapped_df = mapped_df[mapped_df['organ'].notnull()]
        
        print(f"Successfully mapped {len(mapped_df)} out of {len(indications_df)} indications to organs.")
        print(mapped_df.head(20))  # Display first 20 mappings
        
        # Save the results to CSV files
        mapped_df.to_csv('indication_to_organ_mapping_confident.tsv', index=False, sep='\t')
        all_mappings.to_csv('indication_to_organ_mapping_all.tsv', index=False, sep='\t')
        print("Mappings saved to 'indication_to_organ_mapping_confident.tsv'")
        
        # Summary of organ distribution
        if not mapped_df.empty:
            organ_counts = mapped_df['organ'].value_counts()
            print("\nDistribution of indications across organs:")
            print(organ_counts)
            
            # Indications that couldn't be confidently mapped
            unmapped = indications_df[indications_df['organ'].isnull()]
            print(f"\nNumber of indications that couldn't be confidently mapped: {len(unmapped)}")

if __name__ == "__main__":
    main()