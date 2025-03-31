# Benchmarking

Benchmarking logik:
https://docs.google.com/document/d/1e5D6fSlioPEsjZhHv1Cpmvyu5Mxnei5OkHJSsvbRFSc/edit?tab=t.0

## Requirements
Python >= 3.10.14 with pandas < 3.0, matplotlib, seaborn, plotly, statistics and os.

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
  
  At current stage of our instrument testing we use the follow list (81 embryos):

```
3Bal2_eB 5Sachk3_eB 6069_1_K 7493_K 8388_1_K 8388_3_K Aks1_K Aks_K BOC_K1 BOC_e1 BTR_e3 BTR_e9 Bal1_K Bal1_eB Bal3_e Boc_e4 CHR_K1 CHR_e1 Ch-D_e Chal_K1 Chal_e Dik_e Dik_e1 Fed1_e Fisher1_K Gri2_e HAN_K5 HAN_e5 IlI_K3 IlI_e4 Isa1_e K11 K13 Kaz1_e Kaz3_K Kira1_K2 Kond3_K Kond4_K Kov1_e Kov2_e Kra1_e Kra3_K Kra_e Kul_K Kur3_e Kus1_K MAD_K3 MAD_e3 Mar1_K Pan1_K Pash_K2 Pash_K3 Pash_e2 Pash_e3 Sach1_K_2110 Sach2_K Sav4_e Sav4_eB Say3_K Sheg1_K2 Sheg1_e Shen1_K Shen1_e Shen2_K Shen2_e Shen3_K Shen3_e Shur3_K Ton1_e Vla1_K_0705 Vla1_e Vla9_e Vor1_K XAH_K13 XAH_e13 YAK_e4 Zap_K2 Zap_e2 Zap_e3 embryo6 microchip-c 
```
  
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
  Path to the output.
  As default files will be saved in folder *output* at the directory the script is located.


