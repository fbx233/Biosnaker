import pandas as pd
import re
import argparse
import sys
import logging

def read_file(file_path):
    if file_path.endswith('.csv'):
        return pd.read_csv(file_path, header=None)
    elif file_path.endswith('.xlsx'):
        return pd.read_excel(file_path, header=None)
    elif file_path.endswith('.txt'):
        return pd.read_csv(file_path, delimiter='\t', header=None)
    else:
        raise ValueError(f"Unsupported file format: {file_path}")

def extract_key(name):
    match = re.match(r'^([a-zA-Z0-9]+-[a-zA-Z]+)', name)
    return match.group(1) if match else name

def main():
    parser = argparse.ArgumentParser(description='Generate sample information CSV from gene expression matrix')
    parser.add_argument('--input_files', type=str, required=True, help='Path to the input gene expression matrix file')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output sample information CSV file')
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    try:
        logging.debug(f"Reading input file: {args.input_files}")
        gene_expression_matrix = read_file(args.input_files)

        # Check if the file is read correctly
        if gene_expression_matrix.empty:
            logging.error("Input file is empty or not correctly formatted")
            raise ValueError("Input file is empty or not correctly formatted")
        
        logging.debug("Extracting sample names")
        sample_names = gene_expression_matrix.iloc[0, 1:].tolist()
        logging.debug(f"Sample names: {sample_names}")

        logging.debug("Extracting groups")
        group_keys = [extract_key(sample) for sample in sample_names]
        logging.debug(f"Group keys: {group_keys}")

        # Generate unique group labels
        group_mapping = {}
        group_count = 1
        groups = []
        for key in group_keys:
            if key not in group_mapping:
                group_mapping[key] = f"T{group_count:02d}"
                group_count += 1
            groups.append(group_mapping[key])
        
        logging.debug(f"Groups: {groups}")

        # Create sample information dataframe
        sample_info = pd.DataFrame({
            'sample': sample_names,
            'group': groups
        })

        # Write to CSV file
        sample_info.to_csv(args.output_file, index=False)
        logging.debug(f"Output file written: {args.output_file}")
        
    except Exception as e:
        logging.error(f"An error occurred: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
