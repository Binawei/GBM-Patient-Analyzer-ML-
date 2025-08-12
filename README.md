# GBM Patient Data Analyzer

## Overview
This program analyzes Glioblastoma multiforme (GBM) patient data to find patients who:
1. Appear in the RNA-seq expression data
2. Have values for all specified clinical features

## Requirements
- Python 3.6 or higher
- Standard Python libraries (csv, sys, re, collections)

## Input Files Required
Place these files in the same directory as the program:
- `GBM_RNAseqdata_HTSEQ_FKPM.harmonized.txt` - RNA-seq expression data
- `clinical_patient_GBM.txt` - Patient clinical information
- `clinical_followup_GBM.txt` - Follow-up care information  
- `clinical_drug_GBM.txt` - Drug treatment information

## Compilation
No compilation required - Python is interpreted.

## Usage
```bash
python3 gbm_patient_analyzer.py <input_features_file>
```

### Input File Format
Create a text file with one feature name per line, enclosed in quotes:
```
"primary_therapy_outcome_success"
"person_neoplasm_cancer_status"
"age_at_initial_pathologic_diagnosis"
```

### Example Usage
```bash
# Using the provided example input
python3 gbm_patient_analyzer.py example_input.txt

# Using your own feature list
python3 gbm_patient_analyzer.py my_features.txt
```

## Output
The program creates `patient_analysis_results.txt` containing:
- Patient IDs that meet both conditions
- Values for each requested feature for each patient
- Tab-separated format for easy parsing

### Output Format
```
patient ID    0001    2513    ...
"primary_therapy_outcome_success"    Progressive disease    Stable Disease    ...
"age_at_initial_pathologic_diagnosis"    67    82    ...
...
```

## How It Works
1. **Extract Patient IDs**: Parses RNA-seq file headers to extract patient IDs from TCGA barcodes
2. **Read Clinical Data**: Searches all clinical files for patients with the requested features
3. **Match Patients**: Finds patients present in both RNA-seq data and having all clinical features
4. **Generate Output**: Creates formatted results file

## Error Handling
- Missing files: Program exits with error message
- Invalid input format: Skips malformed lines
- No matching patients: Creates output file indicating no results

## Example Run
```bash
$ python3 gbm_patient_analyzer.py example_input.txt
Looking for patients with features: ['primary_therapy_outcome_success', 'person_neoplasm_cancer_status', ...]
Found 528 patients in RNA-seq data
Found clinical data for 592 patients  
Found 45 patients with all required features
Results written to patient_analysis_results.txt
```

## Troubleshooting
- **File not found**: Ensure all data files are in the same directory
- **Permission denied**: Check file permissions
- **Empty results**: Try with fewer required features to test

## Technical Notes
- Patient ID extraction follows TCGA barcode format: `TCGA-XX-XXXX-...` → `XX-XXXX`
- Handles missing values (empty strings, "NA") by excluding them
- Case-sensitive feature matching
- Memory efficient for large datasets