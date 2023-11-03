import pandas as pd
import numpy as np
import os
import sys
import rpy2.robjects as R
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

path = "/Users/thuda/Desktop/Research/14-CRC/" 
#to-do: post-hoc test, graph
#VDJ Count
Checklist_receptors = ['IGK','IGL'] #["TRA","TRB","TRD","TRG","IGH","IGK","IGL"]
paths = ["1-phs001384_timelapse/","2-COAD/","3-CPTAC2/"]

clinical0 = pd.read_csv(path + paths[0] + "/phs001384_Phenotypes.csv",sep=',')
clinical0['histo_stage'] = pd.to_numeric(clinical0['histo_stage'])
clinical0['tumor_having'] = pd.to_numeric(clinical0['tumor_having'])
clinical0 = clinical0.rename(columns={'SAMPLE_ID':'id','histo_stage':'anova'}) #technically, just how samples are being treated, misnomer
clinical0f = clinical0[['id','anova','tumor_having']].copy()
clinical1 = pd.read_csv(path + paths[1]+ "Jclinical2.csv",sep=',') #has files that aren't dl'd, but doesn't matter
clinical1 = clinical1.rename(columns={'PATIENT_ID':'id'}) 
clinical1['anova'] = 5
clinical1f = clinical1[['id','anova']].copy()
clinical2 = pd.read_csv(path + paths[2] + "/clinical.csv",sep=',') #has files that aren't dl'd, but doesn't matter
clinical2 = clinical2.rename(columns={'case_submitter_id':'id'}) 
clinical2['anova'] = 6
clinical2f = clinical2[['id','anova']].copy()
clinical = pd.concat([clinical0f,clinical1f,clinical2f],ignore_index=True)

pc_sheet0 = pd.read_excel(path + paths[0] + "VDJ_Recoveries_p_table.xlsx")
pc_sheet0 = pc_sheet0.rename(columns={'SAMPLE_ID':'id'}) 
pc_sheet1 = pd.read_excel(path + paths[1] + "VDJ_Recoveries_p_table_rnaseq.xlsx")
pc_sheet1 = pc_sheet1.rename(columns={'case_submitter_id':'id'}) 
pc_sheet2 = pd.read_excel(path + paths[2] + "VDJ_Recoveries_p_table_rnaseq.xlsx")
pc_sheet2 = pc_sheet2.rename(columns={'case_submitter_id':'id'}) 
pc_sheet = pd.concat([pc_sheet0,pc_sheet1,pc_sheet2],ignore_index=True)
pc_sheet = pc_sheet.fillna(method='ffill', axis=0)
mean_counts_df = pd.DataFrame(columns=paths)
conditions = []
meanresults = []
VDJfull = pd.DataFrame()
VDJall = pd.DataFrame()

def run_ANOVA(conditions,pdf):
    VDJpivot0 = pdf.groupby(['id'])[conditions[1]].mean().reset_index(name='characteristic') #conditions[2] - which physchem parameter
    VDJpivot = VDJpivot0.merge(clinical[['id','anova']],on='id', how='left')
    VDJpivot['receptor']=conditions[0]
    VDJpivot = VDJpivot[(VDJpivot['anova']<4) | (VDJpivot['anova']==4)] #select which four compared
    groupsizes = VDJpivot.groupby(["anova"]).mean('characteristic')
    semsizes = VDJpivot.groupby(["anova"])['characteristic'].agg(['sem']) #standard error of mean, sem
    finallist = conditions.copy()
    finallist.extend(groupsizes['characteristic'].tolist())
    finallist.extend(semsizes['sem'].tolist())    
    with localconverter(R.default_converter + pandas2ri.converter):
        R_VDJpivot = R.conversion.py2rpy(VDJpivot)
    R_anova = R.r(r'''
        function(df) {        
            options(warn=-1)
            suppressMessages(library(car))
            df$anova <- as.factor(df$anova)
            levels(df$anova) 
            df$anova <- ordered(df$anova, levels = c(1,2,3,4,5,6))
            res.aov <- aov(characteristic~anova, data = df)
            suppressMessages(pvaluearray <- summary(res.aov)[[1]][["Pr(>F)"]][1])
            suppressMessages(library(multcomp))
            post_test <- glht(res.aov,linfct = mcp(anova = "Tukey"))
            suppressMessages(post_testarray <- summary(post_test)[[10]]$pvalues)
            all <- append(pvaluearray, post_testarray)
            }
    ''')
    results=np.asarray(R_anova(R_VDJpivot)) #pop removes nan residual p-value
    finallist.extend(results)
    return finallist

pcresults = []
data = []
siglevel=1

for receptor in Checklist_receptors:
    df1 = pd.DataFrame(pc_sheet[pc_sheet['Receptor']==receptor])
    condition = [receptor,""]
    for column in pc_sheet:
        if column == "id" or column == "Receptor":
            continue
        condition[1] = column
        data = run_ANOVA(condition,df1)
        pcresults.append(data)

Headingsphyschem1 = ["Receptor","Test", "phys-chem - normal_tl", "pc - tubular_tl","pc - villous_tl", "pc - tumor_tl",
                    "pc - COAD", "pc - CPTAC","ANOVA p-value","1-2","1-3","1-4","2-3","2-4","3-4"]
Headingsphyschem = ["Receptor","Test", "phys-chem - normal_tl", "pc - tubular_tl","pc - villous_tl", "pc - tumor",
                    "SEM - normal", "SEM - tubular","SEM - villous", "SEM - tumor","ANOVA p-value","1-2","1-3","1-4","2-3","2-4","3-4"] 
#"1-2","1-3","1-4","1-5","1-6","2-3","2-4","2-5","2-6","3-4","3-5","3-6","4-5","4-6","5-6"

df_pcresults = pd.DataFrame(pcresults).sort_values(by=8,ascending=True)
df_pcresults.columns = Headingsphyschem
with pd.ExcelWriter("TCR_IG_all_ANOVAresults.xlsx") as writer:
    df_pcresults.to_excel(writer,sheet_name="Physchem (ANOVA)")