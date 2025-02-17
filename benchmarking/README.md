# Benchmarking

Benchmarking logik:
https://docs.google.com/document/d/1e5D6fSlioPEsjZhHv1Cpmvyu5Mxnei5OkHJSsvbRFSc/edit?tab=t.0

## Requirements
Python >= 3.10.14 with pandas <3.0, matplotlib, seaborn, plotly, statistics and os.

# Run
```
python3 universal_benchmarker.py -res predictions.csv -ref all_rearrangeements_in_bp_21.01.csv -se list_of_used_embryos -o OUTPUT_DIR
```

# Options:
  -h, --help
  Show this help message and exit.
  
  -res RESULT_PATH, --result_path RESULT_PATH
  Path to your file for analyses (your CNV-caller results).
                        
  -ref REFERENCE_CNV_PATH, --reference_CNV_path REFERENCE_CNV_PATH
  Path to the reference CNV file. (all_rearrangeements_in_bp_21.01.csv)
 
  
  -se SELECTED_EMBRYOS , --selected_embryos SELECTED_EMBRYOS
  Names of selected embryos.
  
  At current stage of our instrument testing we use the follow list (10 embryos):

```
embryo6 BOC_e1 CHR_K1 CHR_e1 IlI_K3 Pash_e3 XAH_e13 YAK_e4 Kra_e BTR_e3 
```
  
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
  Path to the output.
  As default files will be saved in folder *output* at the directory the script is located.


