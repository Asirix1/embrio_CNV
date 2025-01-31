# Benchmarking

Benchmarking logik:
https://docs.google.com/document/d/1e5D6fSlioPEsjZhHv1Cpmvyu5Mxnei5OkHJSsvbRFSc/edit?tab=t.0

## Requirements
Python >= 3.10.14 with pandas, matplotlib, seaborn, plotly, statistics and os.

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
  
  At current stage of our instrument testing we use the follow list (46 embryos):

```
11425_K Zap_e2 Kul2_K Kond4_K YAK_e4 Boc_K4 Kov2_e HAN_e5 XAH_e13 Fuks1_K MAD_K3 Mat2_K Pash_K2 BOC_K1 Kond5_K Kaz3_K BTR_e3 Kira1_K1 Zap_e3 IlI_K3 Ore1_K BOC_e1 CHR_K1 Vla1_e Pash_e3 CHR_e1 Boc_e4 embryo6 Zap_K2 Vla2_e Kra1_e Mog1_e XAH_K13 Mel6_K Gri2_e Pash_K3 Kul_K Shen3_K Sheg1_e Kur3_e Kra_e HAN_K5 Aks_K Sach2_K MAD_e3 Say3_K 
```
  
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
  Path to the output.
  As default files will be saved in folder *output* at the directory the script is located.


