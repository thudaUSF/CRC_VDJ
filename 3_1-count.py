import pandas as pd
import numpy as np
import os
import sys
import rpy2.robjects as R
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

mainpath = "/Users/thuda/Desktop/Research/14-CRC/" 
Checklist_receptors = ["TRA","TRB","TRD","TRG"]
paths = ["1-phs001384_timelapse/","2-COAD/","3-CPTAC2/"]

clinical0 = pd.read_csv(mainpath + paths[0] + "phs001384_Phenotypes.csv",sep=',')
clinical0['histo_stage'] = pd.to_numeric(clinical0['histo_stage'])
clinical0['tumor_having'] = pd.to_numeric(clinical0['tumor_having'])
clinical0 = clinical0.rename(columns={'SAMPLE_ID':'Case ID'}) #technically, just how samples are being treated, misnomer
clinical0f = clinical0[['Case ID','histo_stage','tumor_having']].copy()

clinical1 = pd.read_csv(mainpath + paths[1] + "Jclinical2.csv",sep=',')
samplesheet1 = pd.read_csv(mainpath + paths[1] + "0-Phenotype_HLA/sample_sheet.tsv",sep='\t') #filter to only dl'd files
clinical1 = clinical1.rename(columns={'PATIENT_ID':'Case ID'}) 
clinical1 = clinical1[clinical1['Case ID'].isin(samplesheet1['Case ID'])]
clinical1['histo_stage'] = 5
clinical1f = clinical1[['Case ID','histo_stage']].copy()

clinical2 = pd.read_csv(mainpath + paths[2] +"clinical.csv")
samplesheet2 = pd.read_table(mainpath + paths[2] + "0-Phenotype_HLA/sample_sheet_CRC_CPTAC2.tsv",sep='\t')
clinical2 = clinical2.rename(columns={'case_submitter_id':'Case ID'}) 
clinical2['histo_stage'] = 6
clinical2f = clinical2[['Case ID','histo_stage']].copy()

clinical = pd.concat([clinical0f,clinical1f,clinical2f])
mean_counts_df = pd.DataFrame(columns=paths)
conditions = []
meanresults = []
VDJfull = pd.DataFrame()
VDJall = pd.DataFrame()

def run_ANOVA(conditions,pdf):
    if conditions[1] == "mean":
        VDJpivot0 = pdf.groupby(['Case ID']).size().reset_index(name='characteristic')
        VDJpivot_has = VDJpivot0.merge(clinical[['Case ID','histo_stage']],on='Case ID', how='left')
        VDJpivot_none = clinical[~clinical['Case ID'].isin(VDJpivot_has['Case ID'])] 
        VDJpivot_none = VDJpivot_none[["Case ID","histo_stage"]] 
        VDJpivot_none['characteristic'] = 0 #Add people who had 0 counts to list
        VDJpivot = pd.concat([VDJpivot_has,VDJpivot_none]).dropna(subset=['histo_stage']).reset_index(drop=True)
        VDJpivot['anova'] = VDJpivot['histo_stage'].copy()
        global VDJall
        VDJpivot['receptor']=conditions[0]
        VDJall = pd.concat([VDJall,VDJpivot])
        VDJpivot = VDJpivot[VDJpivot['anova']<5] #just do within group if doing mean count
    elif conditions[1] == "percent":
        VDJpivot = pdf.merge(VDJpivotall[['Case ID','sum_characteristic']],on='Case ID', how='left')
        #if conditions[0] == "TRD": 
            #VDJpivot.to_csv(conditions[0] + conditions[1] + "fi.csv")
            #sys.exit()
        VDJpivot['characteristic'] = VDJpivot['characteristic']/VDJpivot['sum_characteristic']
        VDJpivot = VDJpivot.dropna(subset='characteristic').reset_index(drop=True)
        VDJpivot = VDJpivot[(VDJpivot['anova']<4) | (VDJpivot['anova']==6)]
        #VDJpivot.to_csv(conditions[0] + conditions[1] + ".csv")
    groupsizes = VDJpivot.groupby(["histo_stage"]).mean('characteristic')
    semsizes = VDJpivot.groupby(["histo_stage"])['characteristic'].agg(['sem']) #standard error of mean, sem
    finallist = conditions.copy()    
    finallist.extend(groupsizes['characteristic'].tolist())
    finallist.extend(semsizes['sem'].tolist())
    print(finallist)
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

def run_KW(conditions,pdf):
    VDJpivot = pdf.merge(VDJpivotall[['Case ID','sum_characteristic']],on='Case ID', how='left')
    #if conditions[0] == "TRD": 
        #VDJpivot.to_csv(conditions[0] + conditions[1] + "fi.csv")
        #sys.exit()
    VDJpivot['characteristic'] = VDJpivot['characteristic']/VDJpivot['sum_characteristic']
    VDJpivot = VDJpivot.dropna(subset='characteristic').reset_index(drop=True)
    #VDJpivot.to_csv(conditions[0] + conditions[1] + ".csv")
    groupsizes = VDJpivot.groupby(["histo_stage"]).mean('characteristic')
    finallist = conditions.copy()
    finallist.extend(groupsizes['characteristic'].tolist())
    with localconverter(R.default_converter + pandas2ri.converter):
        R_VDJpivot = R.conversion.py2rpy(VDJpivot)
    Kruskall = R.r(r'''
        function(df) {
            options(warn=-1)
            df$anova <- as.factor(df$anova)
            df$anova <- ordered(df$anova, levels = c(1,2,3,4,5,6))

            pvalue <- kruskal.test(characteristic~anova, data=df)$p.value
            }
    ''')
    posthoc = R.r(r'''
        function(df) {
            options(warn=-1)
            df$anova <- as.factor(df$anova)
            df$anova <- ordered(df$anova, levels = c(1,2,3,4,5,6))

            suppressMessages(library(FSA))
            dunn <- dunnTest(df$characteristic,df$anova,method="bonferroni")$res
            return(dunn) 
            }
    ''')    
    rresults = Kruskall(R_VDJpivot)
    results=np.asarray(rresults)
    if(results<1):
        rposthoc = posthoc(R_VDJpivot)
    with localconverter(R.default_converter + pandas2ri.converter):
        posthoc = R.conversion.rpy2py(rposthoc)
    posthoclist = posthoc["P.adj"]
    finallist.extend(results)
    finallist.extend(posthoclist)
    return finallist

for receptor in Checklist_receptors:
    VDJfull=pd.DataFrame() #cleans it
    for path in paths:
        if path == "1-phs001384_timelapse/":
            VDJ1 = pd.read_csv(mainpath + path + "1-VDJ/" + receptor + "fm.csv",sep=',')
            VDJ1 = VDJ1.rename(columns={'SAMPLE_ID':'Case ID'}) 
            VDJ1 = VDJ1[['Case ID']].copy()
            #VDJ1 = VDJ1.merge(clinical[['Case ID','histo_stage','tumor_having']],on='Case ID', how='left')
            #VDJfull = VDJfull[VDJfull['tumor_having'] == 1]
        else:
            VDJ1 = pd.read_csv(mainpath + path + "1-VDJ/" + receptor + "fmm.csv",sep=',')
            VDJ1 = VDJ1[VDJ1['Chromosome']!='*']
            VDJ1 = VDJ1[VDJ1.Filename.str.contains('rna_seq')]
            VDJ1 = VDJ1[VDJ1['Sample Type'] == 'Primary Tumor']
            VDJ1 = VDJ1.rename(columns={'case_submitter_id':'Case ID'}) 
            VDJ1 = VDJ1[['Case ID']].copy()
        VDJfull = pd.concat([VDJfull,VDJ1])
    data = run_ANOVA([receptor,"mean"],VDJfull)
    meanresults.append(data)
            
percentresults = []
data = []
siglevel=1

VDJpivotall = VDJall.drop(columns=['histo_stage','anova','receptor'])
VDJpivotall = VDJpivotall.groupby(['Case ID']).sum('characteristic').reset_index()
VDJpivotall = VDJpivotall.rename(columns={'characteristic':'sum_characteristic'})

for receptor in Checklist_receptors:
    condition = [receptor,"percent"]
    VJreceptor = pd.DataFrame(VDJall[VDJall['receptor']==receptor])
    data = run_ANOVA(condition,VJreceptor)
    percentresults.append(data)

Headingsmean = ["Receptor","Test", "Mean_counts - normal", "Mean_counts - tubular","Mean_counts - villous", "Mean_counts - tumor", 
                "SEM - normal", "SEM - tubular","SEM - villous", "SEM - tumor", "ANOVA p-value","1-2","1-3","1-4","2-3","2-4","3-4"]
Headingspercent = ["Receptor","Test", "%counts - normal_tl", "%counts - tubular_tl","%counts - villous_tl", "%counts - tumor_tl", 
                   "SEM - %normal", "SEM - %tubular","SEM - %villous", "SEM - %tumor",
                   "ANOVA p-value","1-2","1-3","1-4","2-3","2-4","3-4"]

df_countresults = pd.DataFrame(meanresults).sort_values(by=6,ascending=True)
df_countresults.columns = Headingsmean
df_percentresults = pd.DataFrame(percentresults).sort_values(by=8,ascending=True)
df_percentresults.columns = Headingspercent
print("ji")
with pd.ExcelWriter("Results-11-1.xlsx") as writer:
    df_countresults.to_excel(writer,sheet_name="Receptor counts")
    df_percentresults.to_excel(writer,sheet_name="Receptor percents")