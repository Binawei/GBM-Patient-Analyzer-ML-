#!/usr/bin/env python3

import sys
import csv
import re
from collections import defaultdict

def extract_patient_id_from_sample(sample_name):
    """Extract patient ID from TCGA sample barcode"""
    # Sample format: TCGA-28-2513-01A-01R-1850-01
    # Patient ID is the part after TCGA- up to the third hyphen
    match = re.match(r'TCGA-(\d+)-(\d+)', sample_name)
    if match:
        return f"{match.group(1)}-{match.group(2)}"
    return None

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

def read_clinical_data(clinical_files, required_features):
    """Read clinical data and find patients with all required features"""
    patient_data = defaultdict(dict)
    
    for file_path in clinical_files:
        with open(file_path, 'r') as f:
            # Read header line
            header_line = f.readline().strip()
            headers = [h.strip('"') for h in header_line.split(',')]
            
            # Read data lines
            for line in f:
                values = [v.strip('"') for v in line.strip().split(',')]
                
                # Skip first column
                if len(values) > len(headers):
                    values = values[1:]  # Remove first column
                
                # Create row dictionary
                row = dict(zip(headers, values))
                
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

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 gbm_patient_analyzer.py <input_features_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # File paths
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