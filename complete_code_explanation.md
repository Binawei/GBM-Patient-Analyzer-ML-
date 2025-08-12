# Complete Code Explanation: GBM Patient Analyzer

## Overview
This program analyzes Glioblastoma multiforme (GBM) patient data to find patients who appear in RNA-seq expression data AND have values for all specified clinical features.

## Program Structure

### Imports and Dependencies
```python
#!/usr/bin/env python3

import sys
import csv
import re
from collections import defaultdict
```
**Purpose**:
- `sys`: Command line argument handling and program exit
- `csv`: Parse clinical data files (tab/comma separated)
- `re`: Regular expressions for patient ID extraction from TCGA barcodes
- `defaultdict`: Automatic dictionary creation to prevent KeyErrors

---

## Function 1: `extract_patient_id_from_sample()`

```python
def extract_patient_id_from_sample(sample_name):
    """Extract patient ID from TCGA sample barcode"""
    # Sample format: TCGA-28-2513-01A-01R-1850-01
    # Patient ID is the part after TCGA- up to the third hyphen
    match = re.match(r'TCGA-(\d+)-(\d+)', sample_name)
    if match:
        return f"{match.group(1)}-{match.group(2)}"
    return None
```

### Purpose
Converts full TCGA sample barcodes to patient IDs for matching between files.

### Line-by-Line Explanation
1. **Input**: Full sample barcode like `TCGA-28-2513-01A-01R-1850-01`
2. **Regex Pattern**: `r'TCGA-(\d+)-(\d+)'`
   - `TCGA-` - Literal match
   - `(\d+)` - First capture group: one or more digits (tissue source site)
   - `-` - Literal hyphen
   - `(\d+)` - Second capture group: one or more digits (participant ID)
3. **Extract Groups**: `match.group(1)` = "28", `match.group(2)` = "2513"
4. **Return**: Combined as "28-2513"
5. **Error Handling**: Returns `None` if pattern doesn't match

### Example
- Input: `TCGA-28-2513-01A-01R-1850-01`
- Output: `28-2513`

---

## Function 2: `read_rna_seq_patients()`

```python
def read_rna_seq_patients(rna_file):
    """Read RNA-seq file and extract patient IDs"""
    patients = set()
    with open(rna_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for sample in header:
            patient_id = extract_patient_id_from_sample(sample)
            if patient_id:
                patients.add(patient_id)
    return patients
```

### Purpose
Reads RNA-seq file header to get list of all patients with gene expression data.

### Line-by-Line Explanation
1. **Initialize**: `patients = set()` - Use set to avoid duplicates
2. **Open File**: Read RNA-seq expression file
3. **Read Header**: `f.readline().strip().split('\t')` - First line contains sample names
4. **Process Each Sample**: Loop through column headers
5. **Extract Patient ID**: Use previous function to convert barcode to patient ID
6. **Store Valid IDs**: Only add non-None patient IDs to set
7. **Return**: Set of unique patient IDs

### Example Output
```python
{"02-0047", "28-2513", "06-0125", "15-1444", ...}  # 165 patients total
```

---

## Function 3: `read_clinical_data()`

```python
def read_clinical_data(clinical_files, required_features):
    """Read clinical data and find patients with all required features"""
    patient_data = defaultdict(dict)
    
    for file_path in clinical_files:
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Find the patient barcode - it could be in different columns
                patient_barcode = None
                
                # Check common barcode column names
                for col_name in ['bcr_patient_barcode', 'additional_studies']:
                    if col_name in row and row[col_name] and row[col_name].startswith('TCGA-'):
                        patient_barcode = row[col_name]
                        break
                
                # If not found in named columns, search all columns
                if not patient_barcode:
                    for key, value in row.items():
                        if value and isinstance(value, str) and value.startswith('TCGA-'):
                            patient_barcode = value
                            break
                
                if patient_barcode:
                    # Extract patient ID from barcode (e.g., TCGA-02-0001 -> 02-0001)
                    patient_id = patient_barcode.replace('TCGA-', '')
                    
                    # Check each required feature
                    for feature in required_features:
                        if feature in row and row[feature] and row[feature].strip() and row[feature] != 'NA':
                            patient_data[patient_id][feature] = row[feature].strip()
    
    return patient_data
```

### Purpose
Reads all clinical files and extracts patient data for requested features.

### Detailed Breakdown

#### Data Structure Setup
```python
patient_data = defaultdict(dict)
```
- Creates nested dictionary: `patient_data[patient_id][feature] = value`
- `defaultdict(dict)` automatically creates empty dict for new patients

#### File Processing Loop
```python
for file_path in clinical_files:
```
- Processes 3 files: `clinical_patient_GBM.txt`, `clinical_followup_GBM.txt`, `clinical_drug_GBM.txt`
- Features may be scattered across different files

#### CSV Parsing
```python
reader = csv.DictReader(f)
```
- Treats first row as headers
- Each row becomes dictionary: `{'column_name': 'value', ...}`

#### Patient ID Detection (Method 1)
```python
for col_name in ['bcr_patient_barcode', 'additional_studies']:
    if col_name in row and row[col_name] and row[col_name].startswith('TCGA-'):
        patient_barcode = row[col_name]
        break
```
- Checks known column names first
- Validates: column exists, has value, starts with 'TCGA-'

#### Patient ID Detection (Method 2 - Fallback)
```python
if not patient_barcode:
    for key, value in row.items():
        if value and isinstance(value, str) and value.startswith('TCGA-'):
            patient_barcode = value
            break
```
- Searches ALL columns if not found in known locations
- Handles inconsistent data formats

#### Feature Extraction
```python
for feature in required_features:
    if feature in row and row[feature] and row[feature].strip() and row[feature] != 'NA':
        patient_data[patient_id][feature] = row[feature].strip()
```
- **Four-level validation**:
  1. Column exists
  2. Value not empty/None
  3. Value not just whitespace
  4. Value not 'NA' (missing data marker)

### Example Output
```python
{
    "15-1444": {
        "age_at_initial_pathologic_diagnosis": "0",
        "primary_therapy_outcome_success": "YES"
    },
    "28-2513": {
        "age_at_initial_pathologic_diagnosis": "67"
    }
}
```

---

## Function 4: `find_valid_patients()`

```python
def find_valid_patients(rna_patients, clinical_data, required_features):
    """Find patients that exist in RNA-seq data and have all required clinical features"""
    valid_patients = {}
    
    for patient_id in rna_patients:
        if patient_id in clinical_data:
            patient_features = clinical_data[patient_id]
            # Check if patient has all required features
            if all(feature in patient_features for feature in required_features):
                valid_patients[patient_id] = patient_features
    
    return valid_patients
```

### Purpose
Finds intersection of patients: those in RNA-seq data AND having ALL required clinical features.

### Line-by-Line Explanation
1. **Initialize**: Empty dictionary for valid patients
2. **Loop RNA Patients**: Check each patient from RNA-seq data
3. **Check Clinical Data**: `if patient_id in clinical_data` - Patient has some clinical info
4. **Get Features**: Extract patient's clinical features
5. **Validate Completeness**: `all(feature in patient_features for feature in required_features)`
   - Returns `True` only if patient has ALL requested features
6. **Store Valid**: Only patients meeting both conditions are stored

### Logic Flow
```
RNA-seq patients: {A, B, C, D, E}
Clinical patients: {B, C, D, F, G}
Patients with ALL features: {C, D}
Result: {C, D}  (intersection with complete data)
```

---

## Function 5: `write_output()`

```python
def write_output(output_file, valid_patients, required_features):
    """Write results to output file"""
    with open(output_file, 'w') as f:
        if not valid_patients:
            f.write("No patients found with all required features.\n")
            return
        
        # Write header
        patient_ids = sorted(valid_patients.keys())
        f.write("patient ID")
        for pid in patient_ids:
            f.write(f"\t{pid}")
        f.write("\n")
        
        # Write each feature row
        for feature in required_features:
            f.write(f'"{feature}"')
            for pid in patient_ids:
                value = valid_patients[pid].get(feature, "")
                f.write(f"\t{value}")
            f.write("\n")
```

### Purpose
Creates tab-separated output file with patients as columns and features as rows.

### Output Format Creation

#### Handle Empty Results
```python
if not valid_patients:
    f.write("No patients found with all required features.\n")
    return
```

#### Create Header Row
```python
patient_ids = sorted(valid_patients.keys())
f.write("patient ID")
for pid in patient_ids:
    f.write(f"\t{pid}")
f.write("\n")
```
- Sorts patient IDs for consistent output
- Creates: `patient ID    15-1444    28-2513    ...`

#### Create Feature Rows
```python
for feature in required_features:
    f.write(f'"{feature}"')
    for pid in patient_ids:
        value = valid_patients[pid].get(feature, "")
        f.write(f"\t{value}")
    f.write("\n")
```
- Each feature becomes a row
- Uses `get(feature, "")` to handle missing values gracefully
- Creates: `"age_at_initial_pathologic_diagnosis"    0    67    ...`

### Example Output File
```
patient ID	15-1444	28-2513
"age_at_initial_pathologic_diagnosis"	0	67
"primary_therapy_outcome_success"	YES	Progressive Disease
```

---

## Main Function: `main()`

```python
def main():
    if len(sys.argv) != 2:
        print("Usage: python3 gbm_patient_analyzer.py <input_features_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # File paths (assuming they're in the same directory)
    rna_file = "GBM_RNAseqdata_HTSEQ_FKPM.harmonized.txt"
    clinical_files = [
        "clinical_patient_GBM.txt",
        "clinical_followup_GBM.txt", 
        "clinical_drug_GBM.txt"
    ]
    output_file = "patient_analysis_results.txt"
    
    try:
        # Read required features from input file
        required_features = []
        with open(input_file, 'r') as f:
            for line in f:
                feature = line.strip().strip('"')
                if feature:
                    required_features.append(feature)
        
        print(f"Looking for patients with features: {required_features}")
        
        # Read RNA-seq patient IDs
        rna_patients = read_rna_seq_patients(rna_file)
        print(f"Found {len(rna_patients)} patients in RNA-seq data")
        
        # Read clinical data
        clinical_data = read_clinical_data(clinical_files, required_features)
        print(f"Found clinical data for {len(clinical_data)} patients")
        
        # Find valid patients
        valid_patients = find_valid_patients(rna_patients, clinical_data, required_features)
        print(f"Found {len(valid_patients)} patients with all required features")
        
        # Write output
        write_output(output_file, valid_patients, required_features)
        print(f"Results written to {output_file}")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
```

### Purpose
Orchestrates the entire analysis pipeline and handles user interaction.

### Execution Flow

#### 1. Command Line Validation
```python
if len(sys.argv) != 2:
    print("Usage: python3 gbm_patient_analyzer.py <input_features_file>")
    sys.exit(1)
```
- Ensures exactly one argument (input file) provided
- Shows usage message if incorrect

#### 2. File Path Setup
```python
input_file = sys.argv[1]
rna_file = "GBM_RNAseqdata_HTSEQ_FKPM.harmonized.txt"
clinical_files = ["clinical_patient_GBM.txt", "clinical_followup_GBM.txt", "clinical_drug_GBM.txt"]
output_file = "patient_analysis_results.txt"
```
- Gets input file from command line
- Defines all required data files
- Sets output file name

#### 3. Read Input Features
```python
required_features = []
with open(input_file, 'r') as f:
    for line in f:
        feature = line.strip().strip('"')
        if feature:
            required_features.append(feature)
```
- Reads feature names from input file
- Removes quotes and whitespace
- Skips empty lines

#### 4. Execute Analysis Pipeline
```python
# Step 1: Get RNA-seq patients
rna_patients = read_rna_seq_patients(rna_file)

# Step 2: Get clinical data
clinical_data = read_clinical_data(clinical_files, required_features)

# Step 3: Find intersection
valid_patients = find_valid_patients(rna_patients, clinical_data, required_features)

# Step 4: Generate output
write_output(output_file, valid_patients, required_features)
```

#### 5. Progress Reporting
```python
print(f"Looking for patients with features: {required_features}")
print(f"Found {len(rna_patients)} patients in RNA-seq data")
print(f"Found clinical data for {len(clinical_data)} patients")
print(f"Found {len(valid_patients)} patients with all required features")
print(f"Results written to {output_file}")
```
- Provides user feedback at each step
- Shows counts for debugging

#### 6. Error Handling
```python
except FileNotFoundError as e:
    print(f"Error: File not found - {e}")
    sys.exit(1)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
```
- Handles missing files gracefully
- Catches and reports other errors
- Exits with error code for scripting

---

## Program Entry Point

```python
if __name__ == "__main__":
    main()
```
- Ensures `main()` only runs when script is executed directly
- Allows importing functions without running the program

---

## Complete Execution Example

### Input File (`simple_test.txt`)
```
"age_at_initial_pathologic_diagnosis"
"primary_therapy_outcome_success"
```

### Execution
```bash
python3 gbm_patient_analyzer.py simple_test.txt
```

### Step-by-Step Processing
1. **Read Features**: `["age_at_initial_pathologic_diagnosis", "primary_therapy_outcome_success"]`
2. **RNA-seq Processing**: Extract 165 patient IDs from gene expression file
3. **Clinical Processing**: Find patients with both age and therapy outcome data
4. **Intersection**: Find 1 patient (`15-1444`) with both RNA-seq and complete clinical data
5. **Output**: Create results file with patient data

### Console Output
```
Looking for patients with features: ['age_at_initial_pathologic_diagnosis', 'primary_therapy_outcome_success']
Found 165 patients in RNA-seq data
Found clinical data for 1185 patients
Found 1 patients with all required features
Results written to patient_analysis_results.txt
```

### Output File
```
patient ID	15-1444
"age_at_initial_pathologic_diagnosis"	0
"primary_therapy_outcome_success"	YES
```

---

## Key Design Principles

### 1. Robustness
- Handles inconsistent data formats
- Multiple validation levels
- Graceful error handling

### 2. Efficiency
- Only processes requested features
- Uses sets for fast lookups
- Memory-efficient data structures

### 3. Flexibility
- Works with any number of features
- Handles multiple clinical files
- Extensible architecture

### 4. Data Quality
- Filters out empty/invalid values
- Handles missing data markers
- Validates TCGA barcode formats

### 5. User Experience
- Clear progress reporting
- Helpful error messages
- Standard command-line interface

This program successfully solves the bioinformatics challenge of finding patients with both genomic and complete clinical data, handling the messy realities of real-world medical datasets.