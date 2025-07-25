import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import statistics
import plotly.graph_objects as go
import argparse
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=pd.errors.SettingWithCopyWarning)


def main(preresult_path, reference_CNV_path, selected_embryos, output_dir):


    if os.path.isdir(output_dir):
        print(f"Files will be saved to {output_dir}.")
    else:
        os.mkdir(output_dir)
        print(f"Files will be saved to {output_dir}.")


    preresult=pd.read_csv(preresult_path, sep='\t')
    preresult.rename(columns={'Column': 'Sample (E-embryo, K-biopsy)', 'Chromosome': 'SV_chrom', 'Start': 'SV_start', 'End':'SV_end', 'Class':'Item_count'}, inplace=True)
    preresult.insert(preresult.columns.get_loc('SV_end')+1, 'SV_length', [(preresult['SV_end'][i]-preresult['SV_start'][i]) for i in range(len(preresult))])

    preresult['SV_chrom'].update(str(str(i).removeprefix('chr')) for i in preresult['SV_chrom'])

    preresult.insert(preresult.columns.get_loc('Parameter')+1, 'Status', [0 for _ in range(len(preresult))])


    preresult['Parameter'].update(float(i) for i in preresult['Parameter'])

    quality_list=list(set([i for i in preresult['Parameter']]))
    quality_list.sort()

    result=preresult.head(0)
    for i in range(len(preresult)):
        if i+1!=len(preresult) and preresult['Status'][i]==0 and preresult['Parameter'][i]==preresult['Parameter'][i+1] and preresult['Sample (E-embryo, K-biopsy)'][i]==preresult['Sample (E-embryo, K-biopsy)'][i+1] and preresult['SV_chrom'][i]==preresult['SV_chrom'][i+1] and preresult['Item_count'][i]==preresult['Item_count'][i+1] and (preresult['SV_end'][i]-preresult['SV_start'][i+1])<max(preresult['SV_length'][i],preresult['SV_length'][i+1]):
            preresult['Status'][i+1]=1
            new_row=preresult.loc[[i]]
            new_row.reset_index(drop=True, inplace=True)
            new_row['SV_end'][0]=preresult['SV_end'][i+1]
            while True:
                i+=1
                using=0
                if preresult['Sample (E-embryo, K-biopsy)'][i]==preresult['Sample (E-embryo, K-biopsy)'][i+1] and preresult['Parameter'][i]==preresult['Parameter'][i+1] and preresult['SV_chrom'][i]==preresult['SV_chrom'][i+1] and preresult['Item_count'][i]==preresult['Item_count'][i+1] and (preresult['SV_end'][i]-preresult['SV_start'][i+1])<max(preresult['SV_length'][i],preresult['SV_length'][i+1]):
                    new_row['SV_end'][0]=preresult['SV_end'][i+1]
                    preresult['Status'][i+1]=1
                else:
                    new_row['SV_length'][0]=new_row['SV_end'][0]-new_row['SV_start'][0]
                    result=pd.concat([result, new_row], ignore_index=True)
                    break
        elif preresult['Status'][i]==0:
            new_row=preresult.loc[[i]]
            new_row.reset_index(drop=True, inplace=True)
            result=pd.concat([result, new_row], ignore_index=True)



    control_CNV_pre=pd.read_csv(reference_CNV_path)
    #path to reference rearrangements


    pre_contr_embr_list=list(set(list(control_CNV_pre.query('`All_rearrangements`!="Euploid embryo" and `All_rearrangements`!="without PGT"')['Sample (E-embryo, K-biopsy)'])))


    for j in pre_contr_embr_list:
        vision_table=control_CNV_pre.query(f'`Sample (E-embryo, K-biopsy)`=="{j}"')
        s=0
        for i in list(vision_table['Chromosome']):
            if i not in ['X', 'Y']:
                s=1
        vision_table.reset_index(drop=True, inplace=True)
        if s==0:
            addition=vision_table.head(1)
            addition['All_rearrangements'][0]="Euploid embryo"
            addition['Region_by_ISCN'][0]=pd.NaT
            addition['Region_in_bp'][0]=pd.NaT
            addition['Chromosome'][0]=pd.NaT
            addition['Length'][0]=pd.NaT
            addition['Item_count'][0]=pd.NaT
            addition['Mosaicism_main'][0]=pd.NaT
            addition['Mosaicism_demo'][0]=pd.NaT
            addition['Normal_karyotype'][0]=1
            control_CNV_pre=control_CNV_pre[control_CNV_pre['Sample (E-embryo, K-biopsy)'] != str(j)]
            control_CNV_pre=pd.concat([control_CNV_pre, addition], ignore_index=True)
            control_CNV_pre.reset_index(drop=True, inplace=True)

    #with_mosaicism but then embryos with only mosaic rearrangements are "euplodisied" and other mosaic rearrangements are not accounted for even they are detected or not (no FP and no TP)
    control_CNV=control_CNV_pre.query('`Length`>=10**7 or `All_rearrangements`=="Euploid embryo"').query('`All_rearrangements`=="Euploid embryo" or `Mosaicism_main`<=1 or `Mosaicism_demo`=="2-3" or `Mosaicism_demo`=="1-2"').query('`Chromosome`!="X" and `Chromosome`!="Y"')#.query('`Sequenced Read Pairs`>10000000')
    control_CNV=control_CNV.reset_index(drop=True)
 
    control_CNV.insert(control_CNV.columns.get_loc('Chromosome')+1, 'CNV_start', [str(i).split(':')[-1].split('-')[0] for i in control_CNV['Region_in_bp']])
    control_CNV.insert(control_CNV.columns.get_loc('Chromosome')+2, 'CNV_end', [str(i).split(':')[-1].split('-')[-1] for i in control_CNV['Region_in_bp']])



    embryos_in_comp=list(set([i for i in control_CNV['Sample (E-embryo, K-biopsy)']]))

    used_embryos_in_comp=[]


    for i in (embryos_in_comp):
        if i in list(selected_embryos):
            used_embryos_in_comp.append(i)
    #embryos, wich will be analysed




    #converting Mosaicism_demo (0) and "Euploid embryo" (2) to Mosaicism_main for easier analyse
    for i in range(len(control_CNV)):
        if  control_CNV['Mosaicism_demo'][i]=="2-3" or control_CNV['Mosaicism_demo'][i]=="1-2":
            control_CNV['Mosaicism_main'][i]=0
        if control_CNV['All_rearrangements'][i]=='Euploid embryo':
            control_CNV['Mosaicism_main'][i]=2



    #Euplodysation of embryos with only mosaic rearrangenments

    control_CNV_euplodised=control_CNV.copy()

    for i in range(len(control_CNV)):
        if control_CNV['Mosaicism_main'][i]<1:
            id=1
            for j in control_CNV.query(f'`Sample (E-embryo, K-biopsy)`=="{control_CNV["Sample (E-embryo, K-biopsy)"][i]}"').reset_index(drop=True)['Mosaicism_main']:
                if j==1:
                    #so there are NOT only mosaic rearrangements for this embryo
                    id=0
            if id==1:
                #there are ONLY mosaic rearrangements for this embryo
                control_CNV_euplodised['All_rearrangements'][i]='Euploid embryo'


    q_label=[round(i, 1) for i in quality_list]

    tp_values_list_no_mos=[]
    fp_values_list_no_mos=[]
    fn_values_list_no_mos=[]
    recall_list_no_mos=[]

    precision_list_no_mos=[]

    tp_out_values_list_no_mos=[]
    fp_out_values_list_no_mos=[]
    recall_out_list_no_mos=[]
    precision_out_list_no_mos=[]
    F1_list_no_mos=[]
    F1_out_list_no_mos=[]
    fn_out_values_list_no_mos=[]


    tp_values_list_with_mos=[]
    fp_values_list_with_mos=[]
    fn_values_list_with_mos=[]
    recall_list_with_mos=[]

    precision_list_with_mos=[]

    tp_out_values_list_with_mos=[]
    fp_out_values_list_with_mos=[]
    recall_out_list_with_mos=[]
    precision_out_list_with_mos=[]
    F1_list_with_mos=[]
    F1_out_list_with_mos=[]
    fn_out_values_list_with_mos=[]




    bool_vect=[True]*len(control_CNV)
    for i in range(len(bool_vect)):
        if control_CNV['Sample (E-embryo, K-biopsy)'][i] not in used_embryos_in_comp:
            bool_vect[i]=False
    control_CNV_used=control_CNV.loc[bool_vect].reset_index(drop=True)

    control_CNV_used.query('`All_rearrangements`=="Euploid embryo"').reset_index(drop=True)
    #There are some duplicates. It should be so.
    len(control_CNV_used)
    len(control_CNV_used)
    uniq_CNV=control_CNV_used[control_CNV_used['All_rearrangements']!='Euploid embryo'].reset_index(drop=True)
    uniq_CNV_no_mos=control_CNV_used[control_CNV_used['Mosaicism_main']==1].reset_index(drop=True)
    level_with_mos=len(uniq_CNV)
    #We have 90 CNVs in embryos with not only mosaic CNVs. Total we have 117 mosaik and not mosaic rearrangements (83 not mosaic rearrangements) longer then 10**7 bp and one more shorter (Vla1_e 18q22.3 chr18:71221701-75629292 4.407.591 bp).
    level_no_mos=len(uniq_CNV_no_mos)
    level_no_mos



    used_embryos_in_comp

    metrics_dict_no_mos={}
    
    for quality in quality_list:
        s=0
        table={"embryo": [], 'TP': [], 'FP': [], 'TP': [], 'FN': [], 'Recall': [], 'Precision': [], 'IDs':[], 'Ratio_list':[], 'contr': [], 'preIDs':[]}
        for embryo in used_embryos_in_comp:
            result_main_filt=result.query('`SV_chrom`!="X" and `SV_chrom`!="Y"').query('abs(`SV_length`)>=10**7').query(f"`Sample (E-embryo, K-biopsy)`=='{embryo}'").query(f"`Parameter`=={quality}").reset_index(drop=True)
            #Detected rearrangements for current embryo. Only rearrangements detected on autosoms longer then 10**7 bp are used. query(quality) sets the Threshold level
            embryo_CNV=control_CNV_euplodised.query(f'`Sample (E-embryo, K-biopsy)`=="{embryo}"').reset_index(drop=True)
            #All reference rearrangements for current embryo.
            no_mos_embryo_CNV=embryo_CNV.query('`Mosaicism_main`==1')
            TP=0
            #New TP counter for new embryo.

            P=len(result_main_filt)
            #Positives = detected CNV.
            if len(embryo_CNV)>=0:
                #It haven't any functionality. But let it just stay here. (???If there are any rearrangements? Can it be smaller then 0!??)

                ratio_mean=0
                #How shorter or longer then the reference rearrangemet is the detected rearrangement.
                ratio=[]
                IDs=[]
                #Name of detected rearrangement in Genomenal.
                for contr_n in range(len(embryo_CNV)):
                    TP_rate=0
                    

                    for re_n in range(len(result_main_filt)):

                        if str(embryo_CNV['Region_in_bp'][contr_n])!='nan' and str(embryo_CNV['Region_in_bp'][contr_n])!='NaT':
                            #If it is not realy (without mosaic CNV) Euploid embryo.

                            if embryo_CNV['Mosaicism_main'][contr_n]<1:
                                if (float(result_main_filt['SV_start'][re_n])>=(float(embryo_CNV['CNV_start'][contr_n])-float(embryo_CNV['Length'][contr_n])/2) and float(result_main_filt['SV_start'][re_n])<=(float(embryo_CNV['CNV_start'][contr_n])+float(embryo_CNV['Length'][contr_n])/2)) and (float(result_main_filt['SV_end'][re_n])>=(float(embryo_CNV['CNV_end'][contr_n])-float(embryo_CNV['Length'][contr_n])/2) and float(result_main_filt['SV_end'][re_n])<=(float(embryo_CNV['CNV_end'][contr_n])+float(embryo_CNV['Length'][contr_n])/2)) and int(result_main_filt['SV_chrom'][re_n])==int(embryo_CNV['Chromosome'][contr_n]) and (float(embryo_CNV['Length'][contr_n])*0.5)<abs(float(result_main_filt['SV_length'][re_n])) and int(embryo_CNV['Item_count'][contr_n])==int(result_main_filt['Item_count'][re_n]):
                                    P-=1
                                    
                            elif embryo_CNV['Mosaicism_main'][contr_n]==1:
                                if TP_rate==0 and (float(result_main_filt['SV_start'][re_n])>=(float(embryo_CNV['CNV_start'][contr_n])-float(embryo_CNV['Length'][contr_n])/2) and float(result_main_filt['SV_start'][re_n])<=(float(embryo_CNV['CNV_start'][contr_n])+float(embryo_CNV['Length'][contr_n])/2)) and (float(result_main_filt['SV_end'][re_n])>=(float(embryo_CNV['CNV_end'][contr_n])-float(embryo_CNV['Length'][contr_n])/2) and float(result_main_filt['SV_end'][re_n])<=(float(embryo_CNV['CNV_end'][contr_n])+float(embryo_CNV['Length'][contr_n])/2)) and int(result_main_filt['SV_chrom'][re_n])==int(embryo_CNV['Chromosome'][contr_n]) and (float(embryo_CNV['Length'][contr_n])*0.5)<abs(float(result_main_filt['SV_length'][re_n])) and int(embryo_CNV['Item_count'][contr_n])==int(result_main_filt['Item_count'][re_n]):
                                    TP+=1
                                    TP_rate=1
                                    #there should be only one TP per reference rearrangement

                                    ID=f'{result_main_filt["SV_chrom"][re_n]}_{result_main_filt["SV_length"][re_n]}'
                                    IDs.append(ID)
                                    ratio.append(float(result_main_filt['SV_length'][re_n])/float(embryo_CNV['Length'][contr_n]))
                                    #How shorter or longer then the reference rearrangemet is the detected rearrangement.

                if TP!=0:
                    ratio_mean=statistics.mean(ratio)
                    #How shorter or longer in average then the reference rearrangemets are the detected rearrangements.

            FP=P-TP
            
            FN=len(no_mos_embryo_CNV)-TP
        
            if embryo_CNV['All_rearrangements'][0]=='Euploid embryo':
                recall_t=pd.NaT
                #Recall for each embryo.
                precision_t=pd.NaT
                #Precision for each embryo.
                FN=0
            else:
                recall_t=TP/(TP+FN)
                #Recall for each embryo.
                
                if (TP+FP)==0:
                    precision_t=0
                    #Precision for each embryo if no rearrangemrents are detected.
                    #It's a bad statement, but in this case, I will write so, because it is anyway a technical information for code developer.
                else:
                    #Precision for each embryo if any rearrangemrents are detected.
                    precision_t=TP/(TP+FP)
            




            #Filling the table for current threshold level.
            table['embryo'].append(embryo)
            table['TP'].append(TP)
            table['FP'].append(FP)
            table['FN'].append(FN)
            table['Recall'].append(recall_t)
            table['Precision'].append(precision_t)
            table['Ratio_list'].append(ratio)
            table['IDs'].append(IDs)
            table['contr'].append(embryo_CNV['All_rearrangements'][contr_n])
            table['preIDs'].append(list(result_main_filt['SV_chrom']))
            s+=1

        
        

        metrics=pd.DataFrame(data=table)
	
        metrics=metrics.sort_values(by='embryo')
	
        metrics_dict_no_mos[quality]=metrics.copy()



        if quality==quality_list[0]:
            metrics_start=metrics.copy()
        
        if quality==quality_list[int(len(quality_list)//2)]:
            metrics_middle=metrics.copy()




        fp_values_list_no_mos.append(sum(metrics['FP']))
        #writing FP for current threshold level for ROC-like graph

        tp_values_list_no_mos.append(sum(metrics['TP']))
        #writing TP for current threshold level for ROC-like graph

        fn_values_list_no_mos.append(sum(metrics['FN']))
        #writing FN for current threshold level for ROC-like graph

        recall=sum(metrics['TP'])/(sum(metrics['TP'])+sum(metrics['FN']))
        #calculating CNV-level recall value for current threshhold level

        precision=sum(metrics['TP'])/(sum(metrics['TP'])+sum(metrics['FP']))
        #calculating CNV-level precision value for current threshhold level

        recall_list_no_mos.append(recall)
        #wrinting CNV-level recall value for current threshhold level

        precision_list_no_mos.append(precision)
        #writing CNV-level precision value for current threshhold level


        TP_out=len(metrics.query('`contr`!="Euploid embryo" and (`FP`>0 or `TP`>0)'))
        #calculating embryo-level TP for current threshhold level

        FP_out=len(metrics.query('`contr`=="Euploid embryo" and `FP`>0'))
        #calculating embryo-level FP for current threshhold level

        FN_out=len(metrics.query('`contr`!="Euploid embryo" and (`FP`==0 and `TP`==0)'))
        #calculating embryo-level FN for current threshhold level


        if TP_out!=0:
            recall_out=TP_out/(TP_out+FN_out)
            #calculating CNV-level recall value for current threshhold level
        else:
            #No embryos with rearrangements was detected as TP.    But what if there aren't any rearrangements? Nothing. So, we don't use such datasets.
            recall_out=0

        precision_out=TP_out/(TP_out+FP_out) 
        #calculating embryo-level recall value for current threshhold level
        F1=2*recall*precision/(recall+precision)
        #calculating CNV-level F1 value for current threshhold level

        F1_list_no_mos.append(F1)   
        #writing CNV-level F1 value for current threshhold level

        tp_out_values_list_no_mos.append(TP_out)
        #writing embryo-levels TP for current threshold level for ROC-like graph

        fp_out_values_list_no_mos.append(FP_out)
        #writing embryo-levels FP for current threshold level for ROC-like graph
        
        recall_out_list_no_mos.append(recall_out)
        #wrinting embryo-level recall value for current threshhold level

        precision_out_list_no_mos.append(precision_out)
        #writing CNV-embryo precision value for current threshhold level
        
        F1_out=2*recall_out*precision_out/(recall_out+precision_out)
        #calculating embryo-level F1 value for current threshhold level

        F1_out_list_no_mos.append(F1_out)   
        #writing embryo-level F1 value for current threshhold level

        fn_out_values_list_no_mos.append(FN_out)
        #writing embryo-levels FN for current Threshold level for ROC-like graph

        



    plt.clf()
    y = recall_list_no_mos
    x = q_label
    plt.plot(x, y) 
    plt.ylabel('Recall')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/recall_CNV_no_mos.png')

    plt.clf()
    y = precision_list_no_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('Precision')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/precision_CNV_no_mos.png')
 



    plt.clf()
    y = precision_list_no_mos
    x = recall_list_no_mos
    plt.plot(x, y)
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    plt.savefig(f'{output_dir}/recall_to_precision_CNV_no_mos.png')

    plt.clf()
    y = F1_list_no_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('F1')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/F1_CNV_no_mos.png')
 


    #parameters with best F1 rate
    score={"TP": [tp_values_list_no_mos[F1_list_no_mos.index(max(F1_list_no_mos))]], "FP": [fp_values_list_no_mos[F1_list_no_mos.index(max(F1_list_no_mos))]],  "FN": [fn_values_list_no_mos[F1_list_no_mos.index(max(F1_list_no_mos))]], "Recall": [recall_list_no_mos[F1_list_no_mos.index(max(F1_list_no_mos))]], "Precision": [precision_list_no_mos[F1_list_no_mos.index(max(F1_list_no_mos))]], "F1": [max(F1_list_no_mos)], "Threshold": [q_label[F1_list_no_mos.index(max(F1_list_no_mos))]]}
    sc_t=pd.DataFrame(data=score)
    sc_t.to_csv(f"{output_dir}/no_mos_rate_CNV.csv", index=False)
    sc_t






    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=fp_values_list_no_mos,
        y=tp_values_list_no_mos,
        mode='markers',
        marker=dict(size=10, color='blue'),

        text=[f"TP Samples: {tp_sample}<br>FP Samples: {fp_sample}<br>Threshold: {param}" 
            for tp_sample, fp_sample, param in zip(tp_values_list_no_mos, fp_values_list_no_mos, q_label)],
        hoverinfo='text'
    ))

    fig.update_layout(
        title="TP vs FP",
        xaxis_title="False Positives (FP)",
        yaxis_title="True Positives (TP)",
        showlegend=False,
        width=800,
        height=600,
        font=dict(size=18)
    )

    fig.add_shape(
        type="line",
        x0=0, y0=level_no_mos,  # Start point
        x1=max(fp_values_list_no_mos)+max(fp_values_list_no_mos)/10, y1=level_no_mos,  # End point
        line=dict(color="RoyalBlue", width=2)
    )


    fig.write_html(f"{output_dir}/tp_vs_fp_CNV_no_mos.html")

    plt.clf()

    plt.figure(figsize=(10, 7.5))  # Set the figure size (width, height)

    # Scatter plot
    plt.scatter(
        fp_values_list_no_mos,  # X-axis: False Positives (FP)
        tp_values_list_no_mos,  # Y-axis: True Positives (TP)
        s=100,  # Marker size
        c='blue',  # Marker color
        alpha=0.7  # Transparency
    )

    # Add hover-like annotations (tooltips)
    for fp_sample, tp_sample, param in zip(fp_values_list_no_mos, tp_values_list_no_mos, q_label):
        plt.annotate(
            f"TP: {tp_sample}, FP: {fp_sample}, Thr: {param}",
            xy=(fp_sample, tp_sample),  # Position of the annotation
            xytext=(8, -5),  # Offset for the text
            textcoords='offset points',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.05', fc='yellow', alpha=0.5)  # Background color for the annotation
        )

    # Add a horizontal line
    max_fp = max(fp_out_values_list_no_mos)
    plt.axhline(
        y=level_no_mos,  # Y-value for the line
        color='RoyalBlue',  # Line color
        linestyle='-',  # Line style
        linewidth=2  # Line width
    )

    plt.grid(True, color = "grey", linewidth = "1.4",  axis = "y", linestyle = "-")
    plt.grid(True, color = "grey", linewidth = "1.4",  axis = "x", linestyle = "-")

    # Add labels and title
    plt.title("TP vs FP", fontsize=18)
    plt.xlabel("False Positives (FP)", fontsize=16)
    plt.ylabel("True Positives (TP)", fontsize=16)

    # Adjust font size for tick labels
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # Save the plot as a PNG file
    plt.savefig(f"{output_dir}/tp_vs_fp_CNV_no_mos.png", bbox_inches='tight', dpi=300)



    plt.clf()
    y = recall_out_list_no_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('Recall (Embryo level)')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/recall_embryo_no_mos.png')
 

    plt.clf()
    y = precision_out_list_no_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('Precision (Embryo level)')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/recall_embryo_no_mos.png')
 

    plt.clf()
    y = F1_out_list_no_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('F1')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/F1_embryo_no_mos.png')
 

    plt.clf()
    y = precision_out_list_no_mos
    x = recall_out_list_no_mos
    plt.plot(x, y)
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    plt.savefig(f'{output_dir}/recall_to_precision_embryo_no_mos.png')


    #parameters with best F1 rate
    score_out={"TP": [tp_out_values_list_no_mos[F1_out_list_no_mos.index(max(F1_out_list_no_mos))]], "FP": [fp_out_values_list_no_mos[F1_out_list_no_mos.index(max(F1_out_list_no_mos))]], "FN": [fn_out_values_list_no_mos[F1_out_list_no_mos.index(max(F1_out_list_no_mos))]], "Recall": [recall_out_list_no_mos[F1_out_list_no_mos.index(max(F1_out_list_no_mos))]], "Precision": [precision_out_list_no_mos[F1_out_list_no_mos.index(max(F1_out_list_no_mos))]], "F1": [max(F1_out_list_no_mos)], "Threshold": [q_label[F1_out_list_no_mos.index(max(F1_out_list_no_mos))]]}
    sc_t_out=pd.DataFrame(data=score_out)
    sc_t_out.to_csv(f'{output_dir}/no_mos_rate_embryo.csv', index=False)
    sc_t_out


    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=fp_out_values_list_no_mos,
        y=tp_out_values_list_no_mos,
        mode='markers',
        marker=dict(size=10, color='blue'),

        text=[f"TP Samples: {tp_sample}<br>FP Samples: {fp_sample}<br>Threshold: {param}" 
            for tp_sample, fp_sample, param in zip(tp_out_values_list_no_mos, fp_out_values_list_no_mos, q_label)],
        hoverinfo='text'
    ))

    fig.update_layout(
        title="TP vs FP (Embryo level)",
        xaxis_title="False Positives (FP)",
        yaxis_title="True Positives (TP)",
        showlegend=False,
        width=800,
        height=600,
        font=dict(size=18)
    )

    fig.add_shape(
        type="line",
        x0=0, y0=len(metrics.query('`contr`!="Euploid embryo"')),  # Start point
        x1=max(fp_out_values_list_no_mos)+max(fp_out_values_list_no_mos)/10, y1=len(metrics.query('`contr`!="Euploid embryo"')),  # End point
        line=dict(color="RoyalBlue", width=2)
    )


    fig.write_html(f"{output_dir}/tp_vs_fp_embryo_no_mos.html")


    plt.clf()
    # Create the plot
    plt.figure(figsize=(10, 7.5))  # Set the figure size (width, height)

    # Scatter plot
    plt.scatter(
        fp_out_values_list_no_mos,  # X-axis: False Positives (FP)
        tp_out_values_list_no_mos,  # Y-axis: True Positives (TP)
        s=100,  # Marker size
        c='blue',  # Marker color
        alpha=0.7  # Transparency
    )

    # Add hover-like annotations (tooltips)
    for fp_sample, tp_sample, param in zip(fp_out_values_list_no_mos, tp_out_values_list_no_mos, q_label):
        plt.annotate(
            f"TP: {tp_sample}, FP: {fp_sample}, Thr: {param}",
            xy=(fp_sample, tp_sample),  # Position of the annotation
            xytext=(8, -5),  # Offset for the text
            textcoords='offset points',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.05', fc='yellow', alpha=0.5)  # Background color for the annotation
        )

    # Add a horizontal line
    max_fp = max(fp_out_values_list_no_mos)
    plt.axhline(
        y=len(metrics.query('`contr`!="Euploid embryo"')),  # Y-value for the line
        color='RoyalBlue',  # Line color
        linestyle='-',  # Line style
        linewidth=2  # Line width
    )

    plt.grid(True, color = "grey", linewidth = "1.4",  axis = "y", linestyle = "-")
    plt.grid(True, color = "grey", linewidth = "1.4",  axis = "x", linestyle = "-")

    # Add labels and title
    plt.title("TP vs FP (Embryo level)", fontsize=18)
    plt.xlabel("False Positives (FP)", fontsize=16)
    plt.ylabel("True Positives (TP)", fontsize=16)

    # Adjust font size for tick labels
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # Save the plot as a PNG file
    plt.savefig(f"{output_dir}/tp_vs_fp_embryo_no_mos.png", bbox_inches='tight', dpi=300)



    metrics_start.to_csv(f'{output_dir}/no_mos_metrics_low_threshold.csv', index=False)
    metrics_start
    
    metrics_dict_no_mos[q_label[F1_list_no_mos.index(max(F1_list_no_mos))]].to_csv(f'{output_dir}/no_mos_metrics_optimal_threshold.csv', index=False)



    metrics_middle.to_csv(f'{output_dir}/no_mos_metrics_middle_threshold.csv', index=False)
    metrics_middle



    metrics_high=metrics.copy()
    metrics_high.to_csv(f'{output_dir}/no_mos_metrics_high_threshold.csv', index=False)
    metrics_high

















    bool_vect=[True]*len(control_CNV)
    for i in range(len(bool_vect)):
        if control_CNV['Sample (E-embryo, K-biopsy)'][i] not in used_embryos_in_comp:
            bool_vect[i]=False
    control_CNV_used=control_CNV.loc[bool_vect].reset_index(drop=True)


    control_CNV_used.query('`All_rearrangements`=="Euploid embryo"').reset_index(drop=True)
    #There are some duplicates. It should be so.
    len(control_CNV_used)
    len(control_CNV_used)
    uniq_CNV=control_CNV_used[control_CNV_used['All_rearrangements']!='Euploid embryo'].reset_index(drop=True)
    uniq_CNV_with_mos=control_CNV_used[control_CNV_used['Mosaicism_main']==1].reset_index(drop=True)
    level_with_mos=len(uniq_CNV)

    #We have 90 CNVs in embryos with not only mosaic CNVs. Total we have 117 mosaik and not mosaic rearrangements (83 not mosaic rearrangements) longer then 10**7 bp and one more shorter (Vla1_e 18q22.3 chr18:71221701-75629292 4.407.591 bp).






    used_embryos_in_comp

    metrics_dict_with_mos={}

    for quality in quality_list:
        s=0
        table={"embryo": [], 'TP': [], 'FP': [], 'TP': [], 'FN': [], 'Recall': [], 'Precision': [], 'IDs':[], 'Ratio_list':[], 'contr': [], 'preIDs':[]}
        for embryo in used_embryos_in_comp:
            result_main_filt=result.query('`SV_chrom`!="X" and `SV_chrom`!="Y"').query('abs(`SV_length`)>=10**7').query(f"`Sample (E-embryo, K-biopsy)`=='{embryo}'").query(f"`Parameter`=={quality}").reset_index(drop=True)
            #Detected rearrangements for current embryo. Only rearrangements detected on autosoms longer then 10**7 bp are used. query(quality) sets the Threshold level
            embryo_CNV=control_CNV.query(f'`Sample (E-embryo, K-biopsy)`=="{embryo}"').reset_index(drop=True)
            #All reference rearrangements for current embryo.
            TP=0
            #New TP counter for new embryo.

            P=len(result_main_filt)
            #Positives = detected CNV.
            if len(embryo_CNV)>=0:
                #It haven't any functionality. But let it just stay here. (???If there are any rearrangements? Can it be smaller then 0!??)

                ratio_mean=0
                #How shorter or longer then the reference rearrangemet is the detected rearrangement.
                ratio=[]
                IDs=[]
                #Name of detected rearrangement in Genomenal.
                for contr_n in range(len(embryo_CNV)):
                    TP_rate=0
                    

                    for re_n in range(len(result_main_filt)):

                        if str(embryo_CNV['Region_in_bp'][contr_n])!='nan' and str(embryo_CNV['Region_in_bp'][contr_n])!='NaT':
                            #If it is not realy (without mosaic CNV) Euploid embryo.

                            if TP_rate==0 and (float(result_main_filt['SV_start'][re_n])>=(float(embryo_CNV['CNV_start'][contr_n])-float(embryo_CNV['Length'][contr_n])/2) and float(result_main_filt['SV_start'][re_n])<=(float(embryo_CNV['CNV_start'][contr_n])+float(embryo_CNV['Length'][contr_n])/2)) and (float(result_main_filt['SV_end'][re_n])>=(float(embryo_CNV['CNV_end'][contr_n])-float(embryo_CNV['Length'][contr_n])/2) and float(result_main_filt['SV_end'][re_n])<=(float(embryo_CNV['CNV_end'][contr_n])+float(embryo_CNV['Length'][contr_n])/2)) and int(result_main_filt['SV_chrom'][re_n])==int(embryo_CNV['Chromosome'][contr_n]) and (float(embryo_CNV['Length'][contr_n])*0.5)<abs(float(result_main_filt['SV_length'][re_n])) and int(embryo_CNV['Item_count'][contr_n])==int(result_main_filt['Item_count'][re_n]):
                                TP+=1
                                TP_rate=1
                                #there should be only one TP per reference rearrangement

                                ID=f'{result_main_filt["SV_chrom"][re_n]}_{result_main_filt["SV_length"][re_n]}'
                                IDs.append(ID)
                                ratio.append(float(result_main_filt['SV_length'][re_n])/float(embryo_CNV['Length'][contr_n]))
                                #How shorter or longer then the reference rearrangemet is the detected rearrangement.

                if TP!=0:
                    ratio_mean=statistics.mean(ratio)
                    #How shorter or longer in average then the reference rearrangemets are the detected rearrangements.

            FP=P-TP
            
            FN=len(embryo_CNV)-TP
        
            if embryo_CNV['All_rearrangements'][0]=='Euploid embryo':
                recall_t=pd.NaT
                #Recall for each embryo.
                precision_t=pd.NaT
                #Precision for each embryo.
                FN=0
            else:
                recall_t=TP/(TP+FN)
                #Recall for each embryo.
                
                if (TP+FP)==0:
                    precision_t=0
                    #Precision for each embryo if no rearrangemrents are detected.
                    #It's a bad statement, but in this case, I will write so, because it is anyway a technical information for code developer.
                else:
                    #Precision for each embryo if any rearrangemrents are detected.
                    precision_t=TP/(TP+FP)
            




            #Filling the table for current threshold level.
            table['embryo'].append(embryo)
            table['TP'].append(TP)
            table['FP'].append(FP)
            table['FN'].append(FN)
            table['Recall'].append(recall_t)
            table['Precision'].append(precision_t)
            table['Ratio_list'].append(ratio)
            table['IDs'].append(IDs)
            table['contr'].append(embryo_CNV['All_rearrangements'][contr_n])
            table['preIDs'].append(list(result_main_filt['SV_chrom']))
            s+=1

        


        metrics=pd.DataFrame(data=table)

        metrics=metrics.sort_values(by='embryo')
        
        metrics_dict_with_mos[quality]=metrics.copy()

        if quality==quality_list[0]:
            metrics_start=metrics.copy()
        
        if quality==quality_list[int(len(quality_list)//2)]:
            metrics_middle=metrics.copy()




        fp_values_list_with_mos.append(sum(metrics['FP']))
        #writing FP for current threshold level for ROC-like graph

        tp_values_list_with_mos.append(sum(metrics['TP']))
        #writing TP for current threshold level for ROC-like graph

        fn_values_list_with_mos.append(sum(metrics['FN']))
        #writing FN for current threshold level for ROC-like graph

        recall=sum(metrics['TP'])/(sum(metrics['TP'])+sum(metrics['FN']))
        #calculating CNV-level recall value for current threshhold level

        precision=sum(metrics['TP'])/(sum(metrics['TP'])+sum(metrics['FP']))
        #calculating CNV-level precision value for current threshhold level

        recall_list_with_mos.append(recall)
        #wrinting CNV-level recall value for current threshhold level

        precision_list_with_mos.append(precision)
        #writing CNV-level precision value for current threshhold level


        TP_out=len(metrics.query('`contr`!="Euploid embryo" and (`FP`>0 or `TP`>0)'))
        #calculating embryo-level TP for current threshhold level

        FP_out=len(metrics.query('`contr`=="Euploid embryo" and `FP`>0'))
        #calculating embryo-level FP for current threshhold level

        FN_out=len(metrics.query('`contr`!="Euploid embryo" and (`FP`==0 and `TP`==0)'))
        #calculating embryo-level FN for current threshhold level

        


        if TP_out!=0:
            recall_out=TP_out/(TP_out+FN_out)
            #calculating CNV-level recall value for current threshhold level
        else:
            #No embryos with rearrangements was detected as TP.    But what if there aren't any rearrangements? Nothing. So, we don't use such datasets.
            recall_out=0

        precision_out=TP_out/(TP_out+FP_out) 
        #calculating embryo-level recall value for current threshhold level
        F1=2*recall*precision/(recall+precision)
        #calculating CNV-level F1 value for current threshhold level

        F1_list_with_mos.append(F1)   
        #writing CNV-level F1 value for current threshhold level

        tp_out_values_list_with_mos.append(TP_out)
        #writing embryo-levels TP for current threshold level for ROC-like graph

        fp_out_values_list_with_mos.append(FP_out)
        #writing embryo-levels FP for current threshold level for ROC-like graph
        
        recall_out_list_with_mos.append(recall_out)
        #wrinting embryo-level recall value for current threshhold level

        precision_out_list_with_mos.append(precision_out)
        #writing CNV-embryo precision value for current threshhold level
        
        F1_out=2*recall_out*precision_out/(recall_out+precision_out)
        #calculating embryo-level F1 value for current threshhold level

        F1_out_list_with_mos.append(F1_out)   
        #writing embryo-level F1 value for current threshhold level

        fn_out_values_list_with_mos.append(FN_out)
        #writing embryo-levels FN for current threshold level for ROC-like graph





    plt.clf()
    y = recall_list_with_mos
    x = q_label
    plt.plot(x, y) 
    plt.ylabel('Recall')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/recall_CNV_with_mos.png')
 

    plt.clf()
    y = precision_list_with_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('Precision')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/precision_CNV_with_mos.png')
 



    plt.clf()
    y = precision_list_with_mos
    x = recall_list_with_mos
    plt.plot(x, y)
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    plt.savefig(f'{output_dir}/recall_to_precision_CNV_with_mos.png')

    plt.clf()
    y = F1_list_with_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('F1')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/F1_CNV_with_mos.png')
 


    #parameters with best F1 rate
    score={"TP": [tp_values_list_with_mos[F1_list_with_mos.index(max(F1_list_with_mos))]], "FP": [fp_values_list_with_mos[F1_list_with_mos.index(max(F1_list_with_mos))]], "FN": [fn_values_list_with_mos[F1_list_with_mos.index(max(F1_list_with_mos))]], "Recall": [recall_list_with_mos[F1_list_with_mos.index(max(F1_list_with_mos))]], "Precision": [precision_list_with_mos[F1_list_with_mos.index(max(F1_list_with_mos))]], "F1": [max(F1_list_with_mos)], "Threshold": [q_label[F1_list_with_mos.index(max(F1_list_with_mos))]]}
    sc_t=pd.DataFrame(data=score)
    sc_t.to_csv(f"{output_dir}/with_mos_rate_CNV.csv", index=False)
    sc_t






    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=fp_values_list_with_mos,
        y=tp_values_list_with_mos,
        mode='markers',
        marker=dict(size=10, color='blue'),

        text=[f"TP Samples: {tp_sample}<br>FP Samples: {fp_sample}<br>Threshold: {param}" 
            for tp_sample, fp_sample, param in zip(tp_values_list_with_mos, fp_values_list_with_mos, q_label)],
        hoverinfo='text'
    ))

    fig.update_layout(
        title="TP vs FP",
        xaxis_title="False Positives (FP)",
        yaxis_title="True Positives (TP)",
        showlegend=False,
        width=800,
        height=600,
        font=dict(size=18)
    )

    fig.add_shape(
        type="line",
        x0=0, y0=level_with_mos,  # Start point
        x1=max(fp_values_list_with_mos)+max(fp_values_list_with_mos)/10, y1=level_with_mos,  # End point
        line=dict(color="RoyalBlue", width=2)
    )


    fig.write_html(f"{output_dir}/tp_vs_fp_CNV_with_mos.html")

    plt.clf()

    # Create the plot
    plt.figure(figsize=(10, 7.5))  # Set the figure size (width, height)

    # Scatter plot
    plt.scatter(
        fp_values_list_with_mos,  # X-axis: False Positives (FP)
        tp_values_list_with_mos,  # Y-axis: True Positives (TP)
        s=100,  # Marker size
        c='blue',  # Marker color
        alpha=0.7  # Transparency
    )

    # Add hover-like annotations (tooltips)
    for fp_sample, tp_sample, param in zip(fp_values_list_with_mos, tp_values_list_with_mos, q_label):
        plt.annotate(
            f"TP: {tp_sample}, FP: {fp_sample}, Thr: {param}",
            xy=(fp_sample, tp_sample),  # Position of the annotation
            xytext=(8, -5),  # Offset for the text
            textcoords='offset points',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.05', fc='yellow', alpha=0.5)  # Background color for the annotation
        )

    # Add a horizontal line
    max_fp = max(fp_out_values_list_with_mos)
    plt.axhline(
        y=level_with_mos,  # Y-value for the line
        color='RoyalBlue',  # Line color
        linestyle='-',  # Line style
        linewidth=2  # Line width
    )

    plt.grid(True, color = "grey", linewidth = "1.4",  axis = "y", linestyle = "-")
    plt.grid(True, color = "grey", linewidth = "1.4",  axis = "x", linestyle = "-")

    # Add labels and title
    plt.title("TP vs FP", fontsize=18)
    plt.xlabel("False Positives (FP)", fontsize=16)
    plt.ylabel("True Positives (TP)", fontsize=16)

    # Adjust font size for tick labels
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # Save the plot as a PNG file
    plt.savefig(f"{output_dir}/tp_vs_fp_CNV_with_mos.png", bbox_inches='tight', dpi=300)


    plt.clf()
    y = recall_out_list_with_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('Recall (Embryo level)')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/recall_embryo_with_mos.png')
 

    plt.clf()
    y = precision_out_list_with_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('Precision (Embryo level)')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/recall_embryo_with_mos.png')
 

    plt.clf()
    y = F1_out_list_with_mos
    x = q_label
    plt.plot(x, y)
    plt.ylabel('F1')
    plt.xlabel('Threshold')
    plt.savefig(f'{output_dir}/F1_embryo_with_mos.png')
 

    plt.clf()
    y = precision_out_list_with_mos
    x = recall_out_list_with_mos
    plt.plot(x, y)
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    plt.savefig(f'{output_dir}/recall_to_precision_embryo_with_mos.png')


    #parameters with best F1 rate
    score_out={"TP": [tp_out_values_list_with_mos[F1_out_list_with_mos.index(max(F1_out_list_with_mos))]], "FP": [fp_out_values_list_with_mos[F1_out_list_with_mos.index(max(F1_out_list_with_mos))]], "FN": [fn_out_values_list_with_mos[F1_out_list_with_mos.index(max(F1_out_list_with_mos))]], "Recall": [recall_out_list_with_mos[F1_out_list_with_mos.index(max(F1_out_list_with_mos))]], "Precision": [precision_out_list_with_mos[F1_out_list_with_mos.index(max(F1_out_list_with_mos))]], "F1": [max(F1_out_list_with_mos)], "Threshold": [q_label[F1_out_list_with_mos.index(max(F1_out_list_with_mos))]]}
    sc_t_out=pd.DataFrame(data=score_out)
    sc_t_out.to_csv(f"{output_dir}/with_mos_rate_embryo.csv", index=False)
    sc_t_out


    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=fp_out_values_list_with_mos,
        y=tp_out_values_list_with_mos,
        mode='markers',
        marker=dict(size=10, color='blue'),

        text=[f"TP Samples: {tp_sample}<br>FP Samples: {fp_sample}<br>Threshold: {param}" 
            for tp_sample, fp_sample, param in zip(tp_out_values_list_with_mos, fp_out_values_list_with_mos, q_label)],
        hoverinfo='text'
    ))

    fig.update_layout(
        title="TP vs FP (Embryo level)",
        xaxis_title="False Positives (FP)",
        yaxis_title="True Positives (TP)",
        showlegend=False,
        width=800,
        height=600,
        font=dict(size=18)
    )

    fig.add_shape(
        type="line",
        x0=0, y0=len(metrics.query('`contr`!="Euploid embryo"')),  # Start point
        x1=max(fp_out_values_list_with_mos)+max(fp_out_values_list_with_mos)/10, y1=len(metrics.query('`contr`!="Euploid embryo"')),  # End point
        line=dict(color="RoyalBlue", width=2)
    )


    fig.write_html(f"{output_dir}/tp_vs_fp_embryo_with_mos.html")

    plt.clf()

    # Create the plot
    plt.figure(figsize=(10, 7.5))  # Set the figure size (width, height)

    # Scatter plot
    plt.scatter(
        fp_out_values_list_with_mos,  # X-axis: False Positives (FP)
        tp_out_values_list_with_mos,  # Y-axis: True Positives (TP)
        s=100,  # Marker size
        c='blue',  # Marker color
        alpha=0.7  # Transparency
    )

    # Add hover-like annotations (tooltips)
    for fp_sample, tp_sample, param in zip(fp_out_values_list_with_mos, tp_out_values_list_with_mos, q_label):
        plt.annotate(
            f"TP: {tp_sample}, FP: {fp_sample}, Thr: {param}",
            xy=(fp_sample, tp_sample),  # Position of the annotation
            xytext=(8, -5),  # Offset for the text
            textcoords='offset points',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.05', fc='yellow', alpha=0.5)  # Background color for the annotation
        )

    # Add a horizontal line
    max_fp = max(fp_out_values_list_with_mos)
    plt.axhline(
        y=len(metrics.query('`contr`!="Euploid embryo"')),  # Y-value for the line
        color='RoyalBlue',  # Line color
        linestyle='-',  # Line style
        linewidth=2  # Line width
    )

    plt.grid(True, color = "grey", linewidth = "1.4",  axis = "y", linestyle = "-")
    plt.grid(True, color = "grey", linewidth = "1.4",  axis = "x", linestyle = "-")

    # Add labels and title
    plt.title("TP vs FP (Embryo level)", fontsize=18)
    plt.xlabel("False Positives (FP)", fontsize=16)
    plt.ylabel("True Positives (TP)", fontsize=16)

    # Adjust font size for tick labels
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # Save the plot as a PNG file
    plt.savefig(f"{output_dir}/tp_vs_fp_embryo_with_mos.png", bbox_inches='tight', dpi=300)



    metrics_start.to_csv(f'{output_dir}/with_mos_metrics_low_threshold.csv', index=False)
    metrics_start
    

    metrics_dict_with_mos[q_label[F1_list_with_mos.index(max(F1_list_with_mos))]].to_csv(f'{output_dir}/with_mos_metrics_optimal_threshold.csv', index=False)



    metrics_middle.to_csv(f'{output_dir}/with_mos_metrics_middle_threshold.csv', index=False)
    metrics_middle



    metrics_high=metrics.copy()
    metrics_high.to_csv(f'{output_dir}/with_mos_metrics_high_threshold.csv', index=False)
    metrics_high


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmarking of your CNV-caller.")
    parser.add_argument('-res','--result_path', type=str, required=True, help='Path to your file for analyses (your CNV-caller results)')
    parser.add_argument('-ref', '--reference_CNV_path', type=str, required=True, help='Path to the reference CNV file')
    parser.add_argument('-se','--selected_embryos',default=['3Bal2_eB', '5Sachk3_eB', '6069_1_K', '6069_plus', '7493_K', '8388_1_K', '8388_3_K', 'Aks1_K', 'Aks_K', 'BOC_K1', 'BOC_e1', 'BTR_e3', 'BTR_e9', 'Bag3_K', 'Bal1_K', 'Bal1_eB', 'Bal3_e', 'Boc_e4', 'CHR_K1', 'CHR_e1', 'Ezh2_K3110', 'Ezh3_K', 'Fed1_e', 'Fisher1_K', 'Gri2_e', 'HAN_K5', 'HAN_e5', 'IlI_K3', 'Isa1_e', 'Kaz1_e', 'Kaz3_K', 'Kira1_K2', 'Kond3_K', 'Kond4_K', 'Kov1_e', 'Kov2_e', 'Kra1_e', 'Kra3_K', 'Kra_e', 'Kul_K', 'Kur3_K', 'Kur3_e', 'Kus1_K', 'Lam1_K', 'MAD_K3', 'MAD_e3', 'Mar1_K', 'Ore1_e27', 'Pan1_K', 'Pash_K3', 'Pash_e3', 'Pet1_K', 'Sach1_K_2110', 'Sach1_e', 'Sach2_K', 'Sav4_K', 'Sav4_e', 'Sav4_eB', 'Say3_K', 'Sheg1_K2', 'Sheg1_e', 'Shen1_K', 'Shen1_e', 'Shen2_K', 'Shen2_e', 'Shen3_K', 'Shen3_e', 'Shur1_K_2205', 'Shur3_K', 'Tok1_K', 'Ton1_e', 'Vla1_K_0705', 'Vla1_e', 'Vor1_K', 'XAH_K13', 'XAH_e13', 'YAK_e4', 'Zap_K2', 'Zap_e2', 'Zap_e3', 'embryo6'], nargs='+', required=False, help='Names of selected embryos')
    # default value is the list of 81 analysed embryos
    parser.add_argument('-o', '--output_dir', default=f'{os.path.dirname(os.path.abspath(__file__))}/output', type=str, required=False, help='Path to the output')
    args = parser.parse_args()
    main(args.result_path, args.reference_CNV_path, args.selected_embryos, args.output_dir)

