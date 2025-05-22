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
                print(f"Warning: The file {file_name} is missing the ClassKey or SubClassKey columns")
                continue
                
            # Adduct  filtration：the M+H for TG and DG；M+HCOO for all lipids 
            if 'Adduct' in input_df.columns:
                before_filter = input_df.shape[0]
                condition_tg_dg_h = (input_df['ClassKey'].isin(['TG', 'DG'])) & (input_df['Adduct'].str.strip() == 'M+H')
                condition_hcooh = (input_df['Adduct'].str.strip() == 'M+HCOO')
                combined_condition = condition_tg_dg_h | condition_hcooh
                input_df = input_df[~combined_condition]
                after_filter = input_df.shape[0]

            input_df[ontology_column] = input_df.apply(
                lambda row: get_ontology_value(row[class_key_column], row[subclass_key_column], index_df,
                                               ontology_column),axis=1)
            output_file_path = os.path.join(output_folder, file_name)
            input_df.to_excel(output_file_path, index=False)
def get_ontology_value(class_key, subclass_key, index_df, ontology_column):
    matched_rows = index_df[(index_df['ClassKey'] == class_key) & (index_df['SubClassKey'] == subclass_key)]
    if not matched_rows.empty:
        return matched_rows.iloc[0][ontology_column]
    else:
        return None

input_folder = 'Module2-input-lipidsearch'
index_file = 'Module2-1-index of class.xlsx'  #Default Index file for nomenclature conversion
output_folder = 'Module2-output-1'
process_files(input_folder, index_file, output_folder)
