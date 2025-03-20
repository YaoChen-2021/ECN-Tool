#Conversion of Lipid species from MS-DIAL to LipidSearch
import os
import pandas as pd
def process_files(input_folder, index_file, output_folder):
    index_df = pd.read_excel(index_file)
    class_key_column = 'ClassKey'
    subclass_key_column = 'SubClassKey'
    ontology_column = 'Ontology'
    for file_name in os.listdir(input_folder):
        if file_name.endswith('.xlsx') or file_name.endswith('.xls'):
            file_path = os.path.join(input_folder, file_name)
            input_df = pd.read_excel(file_path)
            if class_key_column not in input_df.columns or subclass_key_column not in input_df.columns:
                print(f"Warning: File {file_name} is missing a ClassKey or SubClassKey column.")
                continue

            if 'Adduct' in input_df.columns:
                input_df = input_df[
                    ~((input_df['ClassKey'].isin(['TG', 'DG'])) & (input_df['Adduct'].str.strip() == 'M+H'))]

            # Match MS-DIAL and LipidSearch
            input_df[ontology_column] = input_df.apply(
                lambda row: get_ontology_value(row[class_key_column], row[subclass_key_column], index_df,
                                               ontology_column),
                axis=1)
            output_file_path = os.path.join(output_folder, file_name)
            input_df.to_excel(output_file_path, index=False)
            print(f"Processing complete: {file_name}")

def get_ontology_value(class_key, subclass_key, index_df, ontology_column):
    matched_rows = index_df[(index_df['ClassKey'] == class_key) & (index_df['SubClassKey'] == subclass_key)]
    if not matched_rows.empty:
        return matched_rows.iloc[0][ontology_column]
    else:
        return None

input_folder = 'Module2-input-lipidsearch'
index_file = 'Module2-1-index of class.xlsx'  #Default Index file name
output_folder = 'Module2-output-1'
process_files(input_folder, index_file, output_folder)
