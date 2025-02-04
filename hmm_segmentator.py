# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from pomegranate import HiddenMarkovModel, State, NormalDistribution
import matplotlib.pyplot as plt


def load_bed_regions(file_path):
    """Load BED regions from a file."""
    bed_df = pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 1, 2], names=["Chromosome", "Start", "End"])
    return bed_df


def is_in_filtered_regions(row, filtered_regions):
    """Check if a region overlaps with filtered regions."""
    overlaps = filtered_regions[
        (filtered_regions["Chromosome"] == row["Chromosome"]) &
        (filtered_regions["Start"] <= row["End"]) &
        (filtered_regions["End"] >= row["Start"])
    ]
    return not overlaps.empty


def combine_rows_and_sum(df, size):
    """Combine rows and sum values."""
    grouped = []
    for chrom, group in tqdm(df.groupby('Chromosome'), desc="Combining rows"):
        for i in range(0, len(group), size):
            subset = group.iloc[i:i + size]
            if not subset.empty:
                combined_row = {
                    'Chromosome': chrom,
                    'Start': subset.iloc[0, 1],
                    'End': subset.iloc[-1, 2],
                }
                for col in subset.columns[3:]:
                    combined_row[col] = subset[col].mean()
                grouped.append(combined_row)
    return pd.DataFrame(grouped)


def normalize(values):
    """Normalize values."""
    mean_val = np.mean(values)
    std_val = np.std(values)
    return (values - mean_val) / std_val if std_val > 0 else values


def create_hmm_with_fixed_transitions(disp):
    """Create HMM with fixed transitions using from_matrix."""
    # Define the states and their emission distributions
    state_1 = State(NormalDistribution(0.5, disp), name="1")
    state_2 = State(NormalDistribution(1, disp), name="2")
    state_3 = State(NormalDistribution(1.5, disp), name="3")

    # Define the fixed transition matrix
    transition_matrix = np.array([
        [9.95438464e-01, 4.56153588e-03, 0.00000000e+00],
        [7.10088511e-05, 9.99791875e-01, 1.37115763e-04],
        [0.00000000e+00, 3.37748155e-03, 9.96622518e-01]
    ])

    # Define start and end probabilities (equal for all states)
    start_probabilities = np.array([1/3, 1/3, 1/3])
    end_probabilities = np.array([1/3, 1/3, 1/3])

    # Create the HMM using from_matrix
    model = HiddenMarkovModel.from_matrix(
        transition_matrix,
        [state_1.distribution, state_2.distribution, state_3.distribution],
        start_probabilities,
        end_probabilities
    )

    # Bake the model to finalize it
    model.bake()
    

    return model


def predict_states_with_fixed_hmm(result_df, disp):
    """Predict states using the fixed HMM."""
    predicted_states = {}

    # Create the HMM with fixed transitions
    model = create_hmm_with_fixed_transitions(disp)

    # Predict states for each column (skip the first 3 columns: Chromosome, Start, End)
    for col in result_df.columns[3:]:
        observations = result_df[col].values.flatten()
        predicted_states[col] = np.array(model.predict(observations))

    return predicted_states


def plot_results(result_df, predicted_states):
    """Plot results."""
    unique_chromosomes = result_df['Chromosome'].unique()
    chromosome_boundaries = [result_df[result_df['Chromosome'] == chrom].index[0] for chrom in unique_chromosomes]
    chromosome_boundaries.append(len(result_df))
    midpoints = [(chromosome_boundaries[i] + chromosome_boundaries[i + 1]) // 2 for i in range(len(chromosome_boundaries) - 1)]

    for idx, (col, states) in enumerate(predicted_states.items()):
        normalized_values = result_df[col].values
        plt.figure(figsize=(12, 6))
        plt.scatter(result_df.index, normalized_values, label=f'Normalized Values from {col}', color='blue', alpha=0.1)
        plt.scatter(result_df.index, states, c=states, cmap='viridis', label='Predicted States', alpha=0.7)
        for boundary in chromosome_boundaries[:-1]:
            plt.axvline(x=boundary, color='black', linestyle='--', alpha=0.7)
        plt.xticks(midpoints, unique_chromosomes, rotation=45, ha='center')
        plt.title(f'Normalized Predicted States and Values for {col}')
        plt.xlabel('Chromosome')
        plt.ylabel('Normalized Values / States')
        plt.legend()
        plt.ylim(-0.5, 2.5)
        plt.close()


def get_predictions_df(result_df, predicted_states, disp):
    """Generates a DataFrame with predictions for the current disp threshold."""
    all_data = []
    for col, states in predicted_states.items():
        df = result_df.copy()
        col_base = col.split('_diff')[0]
        df['Class'] = states + 1
        df['Column'] = col_base
        filtered_df = df.loc[(df['Class'] == 1) | (df['Class'] == 3)].copy()
        filtered_df['Group'] = ((filtered_df['Class'] != filtered_df['Class'].shift()) | 
                                (filtered_df['Chromosome'] != filtered_df['Chromosome'].shift())).cumsum()
        aggregated_df = filtered_df.groupby(['Chromosome', 'Group', 'Class'], as_index=False).agg({
            'Start': 'min',
            'End': 'max',
            'Group': 'size'
        })
        aggregated_df = aggregated_df[aggregated_df['Group'] > 25]
        final_df = aggregated_df[['Chromosome', 'Start', 'End', 'Class']].copy()
        final_df.insert(0, 'Column', col_base)
        final_df['Threshold'] = disp
        all_data.append(final_df)
    
    if all_data:
        return pd.concat(all_data, ignore_index=True)
    else:
        return pd.DataFrame()

def main(prediction_coverage, real_coverage, output_file):
    """Basic function for data processing and saving predictions."""
    prediction_coverage = pd.read_csv(prediction_coverage, sep='\t')
    real_coverage = pd.read_csv(real_coverage)

    real_coverage.rename(columns={real_coverage.columns[0]: 'Chromosome', 
                                  real_coverage.columns[1]: 'Start', 
                                  real_coverage.columns[2]: 'End'}, inplace=True)

    columns_to_process = prediction_coverage.columns.difference(['Chromosome'])
    prediction_coverage[columns_to_process] = prediction_coverage[columns_to_process].clip(lower=0)
    prediction_coverage.columns = real_coverage.columns

    merged_df = pd.merge(real_coverage, prediction_coverage, 
                         on=['Chromosome', 'Start', 'End'], 
                         suffixes=('_real', '_pred'))

    real_columns = [col for col in merged_df.columns if col.endswith('_real')]
    pred_columns = [col for col in merged_df.columns if col.endswith('_pred')]

    if len(real_columns) != len(pred_columns):
        raise ValueError("The number of columns in real_coverage and prediction_coverage do not match!")

    for real_col, pred_col in tqdm(zip(real_columns, pred_columns), total=len(real_columns)):
        new_col_name = real_col.replace('_real', '_diff')
        merged_df[new_col_name] = merged_df.apply(
            lambda row: row[real_col] / row[pred_col] if row[pred_col] >= 1 else 0, axis=1
        )

    merged_df = merged_df.drop(columns=real_columns + pred_columns)

    exclude_regions = load_bed_regions("./data/T2T.excluderanges.bed")
    filtered_regions = exclude_regions

    tqdm.pandas(desc="Filtering regions")
    merged_df = merged_df[~merged_df.progress_apply(lambda row: is_in_filtered_regions(row, filtered_regions), axis=1)]

    result_df = combine_rows_and_sum(merged_df, 50)
    result_df_filtered = result_df[~result_df['Chromosome'].isin(['chrX', 'chrY'])]

    all_predictions = []
    
    for disp in tqdm(np.arange(0.2, 2, 0.2), desc="HMM prediction for different thresholds"):
        predicted_states = predict_states_with_fixed_hmm(result_df_filtered, disp)
        plot_results(result_df_filtered, predicted_states)
        predictions_df = get_predictions_df(result_df_filtered, predicted_states, disp)
        if not predictions_df.empty:
            all_predictions.append(predictions_df)
    
    if all_predictions:
        combined_predictions = pd.concat(all_predictions, ignore_index=True)
        combined_predictions.to_csv(output_file, sep='\t', index=False)
    else:
        print("There are no predictions for saving.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process coverage data and save predictions.")
    parser.add_argument("--prediction_coverage", type=str, required=True, help="Path to the prediction coverage file.")
    parser.add_argument("--real_coverage", type=str, required=True, help="Path to the real coverage file.")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output file.")
    args = parser.parse_args()

    main(args.prediction_coverage, args.real_coverage, args.output_file)