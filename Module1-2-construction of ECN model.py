#Construction of ECN models
import os
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

def linear_func(x, a, b):
    return a * x + b
def quadratic_func(x, a, b, c):
    return a * x ** 2 + b * x + c
def filter_x_by_ontology(ontology, x, y):
    if ontology in ['EtherLPA', 'EtherLPC', 'EtherLPE', 'EtherLPI', 'EtherLPG', 'EtherLPS',
                    'LPA', 'LPC', 'LPE', 'LPI', 'LPG', 'LPS']:
        mask = (x >= 10) & (x <= 23)
    elif ontology == 'SM':
        mask = (x >= 26) & (x <= 50)
    elif ontology == 'TG':
        mask = (x >= 30) & (x <= 72)
    elif ontology in ['DG', 'EtherPA', 'EtherPC', 'EtherPE', 'EtherPI', 'EtherPG', 'EtherPS',
                      'PA', 'PC', 'PE', 'PI', 'PG', 'PS']:
        mask = (x >= 25) & (x <= 48)
    else:
        return x, y
    return x[mask], y[mask]

def process_excel(file_path, output_folder):
    data = pd.read_excel(file_path)
    columns = ['Ontology', 'Double bond number', 'Carbon number', 'RT (min)']
    if not all(col in data.columns for col in columns):
        print(f"File {os.path.basename(file_path)} is missing columns and has been skipped.")
        return

    # Data cleansing: ensure that the 'Carbon number' column contains only numbers.
    data['Carbon number'] = pd.to_numeric(data['Carbon number'], errors='coerce')
    data = data.dropna(subset=['Carbon number'])
    data['Carbon number'] = data['Carbon number'].astype(float)
    results = []
    for ontology, group in data.groupby('Ontology'):
        db_groups = group.groupby('Double bond number')
        all_slopes = []
        for db, db_group in db_groups:
            if len(db_group) < 2:
                continue
            x = np.array(db_group['Carbon number'].values, dtype=float)
            y = np.array(db_group['RT (min)'].values, dtype=float)
            x_filtered, y_filtered = filter_x_by_ontology(ontology, x, y)
            if len(x_filtered) < 2:
                continue
            try:
                popt_linear, _ = curve_fit(linear_func, x_filtered, y_filtered)
                all_slopes.append(popt_linear[0])
            except:
                pass
        avg_slope = np.mean(all_slopes) if all_slopes else None
        for db, db_group in db_groups:
            x = np.array(db_group['Carbon number'].values, dtype=float)
            y = np.array(db_group['RT (min)'].values, dtype=float)
            x_filtered, y_filtered = filter_x_by_ontology(ontology, x, y)
            if len(x_filtered) < 2:
                continue
            # Initial linear fitting
            residuals_linear = np.inf
            try:
                popt_linear, _ = curve_fit(linear_func, x_filtered, y_filtered)
                residuals_linear = y_filtered - linear_func(x_filtered, *popt_linear)
                r2_linear = 1 - (np.sum(residuals_linear**2) / np.sum((y_filtered - np.mean(y_filtered))**2))
                if avg_slope and abs(popt_linear[0] - avg_slope) > 0.6:  # Slope constraint
                    continue
            except:
                r2_linear = 0
            # Failed to fit linear equation, used raw data to fit quadratic equations.
            x_filtered_quad, y_filtered_quad = filter_x_by_ontology(ontology, x, y)
            r2_quad = 0
            try:
                if len(x_filtered_quad) >= 3:
                    popt_quad, _ = curve_fit(quadratic_func, x_filtered_quad, y_filtered_quad)
                    residuals_quad = y_filtered_quad - quadratic_func(x_filtered_quad, *popt_quad)
                    r2_quad = 1 - (np.sum(residuals_quad**2) / np.sum((y_filtered_quad - np.mean(y_filtered_quad))**2))
            except:
                pass
            # Remove outliers and refit.
            while (r2_quad < 0.99 and r2_linear < 0.99) and len(x_filtered_quad) > 2:
                residuals = residuals_quad if r2_quad < r2_linear else residuals_linear
                max_residual_idx = np.argmax(np.abs(residuals))
                x_filtered_quad = np.delete(x_filtered_quad, max_residual_idx)
                y_filtered_quad = np.delete(y_filtered_quad, max_residual_idx)
                # Update the linear fit and relax the linear constraints
                try:
                    popt_linear, _ = curve_fit(linear_func, x_filtered_quad, y_filtered_quad)
                    residuals_linear = y_filtered_quad - linear_func(x_filtered_quad, *popt_linear)
                    r2_linear = 1 - (np.sum(residuals_linear**2) / np.sum((y_filtered_quad - np.mean(y_filtered_quad))**2))
                    if avg_slope and abs(popt_linear[0] - avg_slope) > 0.6:
                        r2_linear = 0
                except:
                    r2_linear = 0
                # Update the quadratic equations fitting.
                try:
                    if len(x_filtered_quad) >= 3:
                        popt_quad, _ = curve_fit(quadratic_func, x_filtered_quad, y_filtered_quad)
                        residuals_quad = y_filtered_quad - quadratic_func(x_filtered_quad, *popt_quad)
                        r2_quad = 1 - (np.sum(residuals_quad**2) / np.sum((y_filtered_quad - np.mean(y_filtered_quad))**2))
                    else:
                        r2_quad = 0
                except:
                    pass
            fit_data = {
                'Ontology': ontology,
                'Double bond number': db,
                'X_data': x.tolist(),
                'Y_data': y.tolist(),
                'X_R_data': x_filtered_quad.tolist(),
                'Y_R_data': y_filtered_quad.tolist()}
            if r2_linear >= 0.99 and avg_slope and abs(popt_linear[0] - avg_slope) <= 0.6:
                Equation = f"y = {popt_linear[0]:.4f}x + {popt_linear[1]:.4f}"
                results.append({**fit_data, 'Fit Type': 'Linear', 'Equation': Equation, 'R^2': round(r2_linear, 3)})
            elif r2_quad >= 0.99:
                Equation = f"y = {popt_quad[0]:.4f}x^2 + {popt_quad[1]:.4f}x + {popt_quad[2]:.4f}"
                results.append({**fit_data, 'Fit Type': 'Quadratic', 'Equation': Equation, 'R^2': round(r2_quad, 3)})
    output_file = os.path.join(output_folder, os.path.basename(file_path).replace('.xlsx', '_processed.xlsx'))
    if results:
        results_df = pd.DataFrame(results,
                                  columns=['Ontology', 'Double bond number', 'X_data', 'Y_data', 'X_R_data',
                                           'Y_R_data', 'Fit Type', 'Equation', 'R^2'])
        results_df.to_excel(output_file, sheet_name='Fit Results', index=False)
        print(f"Processed: {os.path.basename(file_path)}, saved: {output_file}")
    else:
        print(f"No fitting results for: {os.path.basename(file_path)}")

def process_folder(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    for file_name in os.listdir(input_folder):
        if file_name.endswith('.xlsx'):
            file_path = os.path.join(input_folder, file_name)
            process_excel(file_path, output_folder)

input_folder = 'Module1-output-1'
output_folder = 'Module1-output-2'
process_folder(input_folder, output_folder)
