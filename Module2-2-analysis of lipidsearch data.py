#Application of ECN model in lipid annotation using LipidSearch
import os
import pandas as pd
import re
import numpy as np
def process_lipid_group_key(df):
    if 'LipidGroupKey' in df.columns:
        carbon_number = []
        double_bond_number = []
        for key in df['LipidGroupKey']:
            if isinstance(key, str):
                match = re.search(r"\D*(\d+):(\d+)\D*", key)
                if match:
                    carbon_number.append(match.group(1))
                    double_bond_number.append(match.group(2))
                else:
                    carbon_number.append(None)
                    double_bond_number.append(None)
            else:
                carbon_number.append(None)
                double_bond_number.append(None)
        df.insert(df.columns.get_loc('LipidGroupKey') + 1, 'Carbon number', carbon_number)
        df.insert(df.columns.get_loc('Carbon number') + 1, 'Double bond number', double_bond_number)
    return df

def solve_Equation(Equation, carbon_number_value):
    try:
        Equation = Equation.replace(" ", "")
        if '=' in Equation:
            Equation = Equation.split('=')[1]
        Equation = Equation.replace("+-", "-")
        Equation = Equation.replace("x", f"*{carbon_number_value}")
        Equation = Equation.replace("^", "**")
        result = eval(Equation)
        return result
    except Exception as e:
        print(f"Error parsing Equation {Equation}: {e}")
        return np.nan

def read_index_table(index_file):
    index_table = pd.read_excel(index_file)
    index_table.columns = index_table.columns.str.strip()
    if 'Ontology' not in index_table.columns:
        print("Error: Index table lacks 'Ontology' column. Skipping further processing.")
        return pd.DataFrame()
    index_table = index_table.dropna(subset=['Ontology'])
    return index_table

def process_excel_step2(input_file, index_table):
    if index_table.empty:
        print("Index table is empty. Skipping file:", input_file)
        return pd.DataFrame()
    try:
        df_input = pd.read_excel(input_file)
        df_input = df_input.dropna(subset=['Ontology'])
        if df_input.empty:
            print(f"File {input_file} contains no valid data after filtering 'Ontology'. Skipping file.")
            return pd.DataFrame()
        df_input = process_lipid_group_key(df_input)
        df_input['Ontology'] = df_input['Ontology'].astype(str)
        df_input['Double bond number'] = pd.to_numeric(df_input['Double bond number'], errors='coerce').fillna(
            0).astype(int)
        df_input['Carbon number'] = pd.to_numeric(df_input['Carbon number'], errors='coerce').fillna(0).astype(int)
        index_table['Ontology'] = index_table['Ontology'].astype(str)
        index_table['Double bond number'] = pd.to_numeric(index_table['Double bond number'], errors='coerce').fillna(
            0).astype(int)
        full_match = pd.merge(df_input, index_table, on=['Ontology', 'Double bond number'], how='left',
                              suffixes=('', '_matched'))
        full_match = full_match.dropna(subset=['Equation'])
        full_match['TheorRT(min)'] = full_match.apply(lambda row: solve_Equation(row['Equation'], row['Carbon number'])
        if pd.notna(row['Equation']) else np.nan, axis=1)
        full_match['δRT(%)'] = ((full_match['TopRT'] - full_match['TheorRT(min)']) / full_match['TopRT']) * 100
        full_match['TheorRT(min)'] = full_match['TheorRT(min)'].round(2)
        full_match['δRT(%)'] = full_match['δRT(%)'].round(2)
        return full_match
    except Exception as e:
        print(f"File {input_file} processing failed with error: {e}")
        return pd.DataFrame()

def filter_and_save(dataframe, output_path):
    try:
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            sheet_data = dataframe
            if "δRT(%)" in sheet_data.columns and "LipidMolec" in sheet_data.columns:
                filtered_data = sheet_data[(sheet_data['δRT(%)'] > -5) & (sheet_data['δRT(%)'] < 5)] #Filter data with δRT(%) >-5 and <5
                processed_data = filtered_data.loc[filtered_data.groupby("LipidMolec")['δRT(%)'].idxmin()].reset_index(drop=True)
                processed_data.to_excel(writer, sheet_name="ProcessedData", index=False)
            else:
                print(f"Required columns 'ΔRT' and 'LipidMolec' not found. Saving original data.")
                sheet_data.to_excel(writer, sheet_name="ProcessedData", index=False)
    except Exception as e:
        print(f"Error saving filtered data to file {output_path}: {e}")

def process_excel_files(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for file_name in os.listdir(input_folder):
        if file_name.endswith('.xlsx'):
            input_path = os.path.join(input_folder, file_name)
            output_path = os.path.join(output_folder, file_name)
            df = pd.read_excel(input_path)
            if 'ClassKey' in df.columns and 'LipidMolec' in df.columns and 'Intensity' in df.columns:
                tg_df = df[df['ClassKey'].str.contains('TG', na=False)].copy()
                def parse_lipid_molec(molec):
                    if pd.isna(molec):
                        return None
                    if '(' in molec and ')' in molec:
                        content = molec.split('(')[-1].split(')')[0]
                        parts = content.split('_')
                        return tuple(sorted(parts))
                    return None
                tg_df['Parsed'] = tg_df['LipidMolec'].apply(parse_lipid_molec)
                filtered_rows = []
                duplicates = []
                for parsed_value, group in tg_df.groupby('Parsed'):
                    if parsed_value is not None:
                        if len(group) > 1:
                            duplicates.append(parsed_value)
                        max_row = group.loc[group['Intensity'].idxmax()]
                        filtered_rows.append(max_row)
                filtered_tg_df = pd.DataFrame(filtered_rows)
                df = pd.concat([df[~df.index.isin(tg_df.index)], filtered_tg_df])
                def process_lipid_class(df, class_key):
                    if not df['ClassKey'].str.contains(class_key, na=False).any():
                        return df
                    class_df = df[df['ClassKey'].str.contains(class_key, na=False)].copy()
                    class_df['Parsed'] = class_df['LipidMolec'].apply(parse_lipid_molec)
                    filtered_class_rows = []
                    class_duplicates = []
                    for parsed_value, group in class_df.groupby('Parsed'):
                        if parsed_value is not None:
                            if len(group) > 1:
                                class_duplicates.append(parsed_value)
                            max_row = group.loc[group['Intensity'].idxmax()]
                            filtered_class_rows.append(max_row)
                    filtered_class_df = pd.DataFrame(filtered_class_rows)
                    return pd.concat([df[~df.index.isin(class_df.index)], filtered_class_df])
                lipid_classes = ['PC', 'PE', 'DG', 'PA', 'PI', 'PG', 'PS']
                for lipid_class in lipid_classes:
                    df = process_lipid_class(df, lipid_class)
            df.to_excel(output_path, index=False)

def main_process(input_folder, index_file, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    index_table = read_index_table(index_file)
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".xlsx"):
            input_path = os.path.join(input_folder, file_name)
            file_name_without_extension = os.path.splitext(file_name)[0]
            output_file_name = f"{file_name_without_extension}_processed.xlsx"
            output_path = os.path.join(output_folder, output_file_name)
            try:
                matched_data = process_excel_step2(input_path, index_table)
                if isinstance(matched_data, pd.DataFrame) and not matched_data.empty:
                    filter_and_save(matched_data, output_path)
                    process_excel_files(output_folder, output_folder)
                else:
                    print(f"{file_name} contains no valid data.")
            except Exception as e:
                print(f"Error processing {file_name}: {e}")

input_folder = "Module2-output-1"
index_file = "Module1-output-2/QE-POS-20MIN-Mixture_processed.xlsx" #The name and path of the ECN model file generated by Module 1-2.
output_folder = "Module2-output-2"
main_process(input_folder, index_file, output_folder)