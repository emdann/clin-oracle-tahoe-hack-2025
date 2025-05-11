# ClinOracle
**Tahoe-DeepDive Hackathon 2025**

## Team Members
- Emma Dann
- Tony Zen
- Ross Giglio
- Kevin Hoffer-Hawlik
- Meer Mustafa

## Project
### Pharmacotranscriptomic representations to predict clinical trial success

### Overview
Large _in vitro_ perturbation screens like Tahoe-100M allow for assessing whether transcriptional responses are predictive of metrics of clinical success like drug approval.

### Motivation
Despite rigorous research efforts, clinical success and drug approval is challenging and difficult to predict. 

### Methods
#### Clinical trial information
We collected clinical trial results associated with the chemical agents screened in Tahoe-100M, annotated for test and approval for the organ of disease.
#### Classifier

#### E-distance


#### LDVAE

#### MrVI
Using the pseudobulked Tahoe-100M data, we trained a MrVI model with sample defined as cell_drug with the union of highly variable genes within cell line as features. We generated two-latent embeddings, the 10-dimensional u-space and the 30-dimensional z-space that were used as input to the classifier.

### Results
#### Code
- `representations/` - scripts to train and transcriptome effect representations on Tahoe-100M 
- `clinical_data_curation/` - scripts to curate clinical trial data
- `approval_prediction_benchmark.ipynb` - benchmark on clinical approval prediction
- `classifier.py` - Benchmarking classifier implementation

#### Data
- `clinical_evidence_data/` - Curated clinical evidence data on Tahoe drugs
- `data_for_classifier/` - input data for benchmarks
- `data/` - misc processed data

#### Discussion and Future Work
