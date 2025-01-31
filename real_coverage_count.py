import argparse
import pyBigWig
import pandas as pd
import numpy as np
from tqdm import tqdm

def load_samples(sample_file):
    """Loads sample names and their corresponding bigWig file paths."""
    samples = {}
    with open(sample_file, 'r') as f:
        for line in f:
            name, bw_path = line.strip().split('\t')
            samples[name] = pyBigWig.open(bw_path)
    return samples

def load_regions(region_file):
    """Loads genomic regions from a CSV file."""
    return pd.read_csv(region_file, header=None, names=["chrom", "start", "end"])

def compute_coverage(samples, regions):
    """Computes total coverage for each sample in each region using numpy."""
    results = []
    for _, row in tqdm(regions.iterrows(), total=len(regions), desc="Computing coverage"):
        chrom, start, end = row
        coverage = [chrom, start, end]
        for sample_name, bw in samples.items():
            try:
                signal = bw.values(chrom, start, end, numpy=True)
                total_coverage = np.nansum(signal) if signal is not None else 0
                coverage.append(total_coverage)
            except Exception:
                coverage.append(0)  # In case of an error, record 0
        results.append(coverage)
    return results

def save_results(output_file, regions, samples, coverage_data):
    """Saves the computed coverage results to an output CSV file."""
    columns = ["chrom", "start", "end"] + list(samples.keys())
    df = pd.DataFrame(coverage_data, columns=columns)
    df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description="Compute bigWig coverage in specified regions.")
    parser.add_argument("--samples", required=True, help="File with sample names and bigWig file paths.")
    parser.add_argument("--regions", required=True, help="File with genomic regions in CSV format.")
    parser.add_argument("--output", required=True, help="Output file name.")
    args = parser.parse_args()
    
    samples = load_samples(args.samples)
    regions = load_regions(args.regions)
    coverage_data = compute_coverage(samples, regions)
    save_results(args.output, regions, samples, coverage_data)
    
    # Close bigWig files
    for bw in samples.values():
        bw.close()

if __name__ == "__main__":
    main()
