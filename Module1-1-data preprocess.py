#Preprocess of MS-DIAL lipid data
import os
import pandas as pd
import re
def process_lipid_group_key(df):
    if 'Name' in df.columns:
        carbon_number = []
        double_bond_number = []
        for key in df['Name']:
            if isinstance(key, str) and ':' in key:
                match = re.search(r"(\d+):(\d+)", key)
                if match:
                    carbon_number.append(match.group(1))
                    double_bond_number.append(match.group(2))
                else:
                    carbon_number.append(key)
                    double_bond_number.append(key)
            else:
                carbon_number.append(key)
                double_bond_number.append(key)
        df.insert(df.columns.get_loc('Name') + 1, 'Carbon number', carbon_number)
        df.insert(df.columns.get_loc('Carbon number') + 1, 'Double bond number', double_bond_number)
    return df

def filter_by_max_height(df):
    required_columns = ['Ontology', 'Carbon number', 'Double bond number', 'Height']
    if not all(col in df.columns for col in required_columns):
        print(f"File is missing required columns, skip processing.")
        return pd.DataFrame()
    grouped = df.groupby('Ontology', group_keys=False)
    processed_data = []
    for _, group in grouped:
        sub_grouped = group.groupby(['Carbon number', 'Double bond number'], group_keys=False)
        for _, sub_group in sub_grouped:
            if len(sub_group) > 1:
                max_row = sub_group[sub_group['Height'] == sub_group['Height'].max()]
                processed_data.append(max_row)
            else:
                processed_data.append(sub_group)
    return pd.concat(processed_data, ignore_index=True)

def main_process(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for filename in os.listdir(input_folder):
        if filename.endswith('.xlsx') or filename.endswith('.xls'):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, filename)
            df = pd.read_excel(input_file)
            df_with_lipid_info = process_lipid_group_key(df)
            processed_df = filter_by_max_height(df_with_lipid_info)
            if not processed_df.empty:
                processed_df.to_excel(output_file, index=False)
                print(f"Processed file {filename} and saved to {output_file}")
            else:
                print(f"{filename} could not be processed successfully, necessary column may be missing.")

input_folder = "Module1-input-MSDIAL"
output_folder = "Module1-output-1"
main_process(input_folder, output_folder)
