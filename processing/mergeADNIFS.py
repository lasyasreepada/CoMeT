import pandas as pd
import numpy as np

""" DATA IO """
# Set directory
dir = '~/Projects/CoMeT/'

# Read in ADNIMERGE
adniMerge = pd.read_csv(dir + 'data/ADNI/participants/ADNIMERGE.csv',dtype={'RID':int, 'VISCODE':str}, parse_dates=['EXAMDATE'],low_memory=False)

# Read in ADNI demographics
adniDems = pd.read_csv(dir + 'data/ADNI/participants/ADNI_Demographics.csv')

# Read in ADNI DNAm
adniDNAm = pd.read_csv(dir + 'data/ADNI/methylation/clocks.csv',parse_dates=['Edate'])

# Read in ADNI Imaging
adniLong = pd.read_csv(dir + 'data/ADNI/imaging/ADNIFSLong.csv',parse_dates=['SCANDATE'],low_memory=False)

# Read in ADNI Cognitive Data
adniCog = pd.read_csv('~/Downloads/ADSP_ADNI_Cognition_Dec2023/ADSP_PHC_COGN_Dec2023.csv',parse_dates=['EXAMDATE'])
adniAVLT = pd.read_csv('~/Downloads/NEUROBAT_10Jun2024.csv',low_memory=False)
adniAVLT = adniAVLT[['RID','VISCODE2','AVDELTOT','AVDELERR2']]

""" METHYLATION DATA """
# Randomly select one replicate per ID per visit
adniDNAmNoReps = adniDNAm.dropna(subset=['RID']).groupby(['RID','Edate']).apply(lambda x: x.sample(1, random_state=42)).reset_index(drop=True)

# Rename Edate to EXAMDATE to enable easy merge 
adniDNAmNoReps.drop(columns=['Unnamed: 0'],inplace=True)
adniDNAmNoReps.rename(columns={'Edate':'EXAMDATE'},inplace=True)

""" COGNITIVE DATA """
dnamCog = adniCog.merge(adniDNAmNoReps, how='outer', on=['RID','EXAMDATE']).rename(columns={'Age':'AGE_SAMPLE'})
dnamCog = dnamCog.merge(adniAVLT,how='outer',on=['RID','VISCODE2'])

""" DEMOGRAPHICS DATA """
adniDems['PTDOBDD'] = 1
adniDems['PTDOB'] = pd.to_datetime(dict(year=adniDems.PTDOBYY, month=adniDems.PTDOBMM, day=adniDems.PTDOBDD))
adniDems = adniDems[['RID','PTDOB']]
dob = adniDems.groupby('RID')['PTDOB'].unique().to_frame().reset_index().explode('PTDOB')
dob = dob.sort_values('PTDOB').groupby('RID').head(1)  

# Extract additional demographic features
gender = adniMerge.groupby('RID')['PTGENDER'].unique().to_frame().reset_index().explode('PTGENDER')
educat = adniMerge.groupby('RID')['PTEDUCAT'].unique().to_frame().reset_index().explode('PTEDUCAT')
apoe = adniMerge.groupby('RID')['APOE4'].unique().to_frame().reset_index().explode('APOE4')

# Create participants dataframe and merge in demographics 
df_participants = pd.DataFrame(adniMerge['RID'].unique(),columns=['RID'])
df_participants = df_participants.merge(gender,on='RID',how='left',validate='1:1')
df_participants = df_participants.merge(educat,on='RID',how='left',validate='1:1')
df_participants = df_participants.merge(apoe,on='RID',how='left',validate='1:1')
df_participants = df_participants.merge(dob,on='RID',how='left',validate='1:1')

# Merge demographic info into longitudinal visits df
df_visits = adniMerge[['RID','VISCODE','EXAMDATE','AGE','DX','FBB','AV45','PIB','ABETA','CDRSB','FLDSTRENG']]
df_visits = df_visits.merge(df_participants,on='RID',how='left',validate='m:1')
df_visits.rename(columns={'AGE':'AGE_BL'},inplace=True)
df_visits['AGE'] = (df_visits['EXAMDATE'] - df_visits['PTDOB']).dt.days / 365.25

""" CLINICAL DATA INTERPOLATION """
# Match and merge Diagnosis
diagnosis = df_visits[['RID','AGE','DX']].sort_values(['RID', 'AGE'])
diagnosis.dropna(subset=['RID','AGE','DX'], how='any', inplace=True)

# Get first and last diagnosis
first_diagnosis = diagnosis.loc[diagnosis.groupby('RID').AGE.idxmin()].rename(columns={'AGE':'Age_DX_BL', 'DX':'DX_BL'})
last_diagnosis = diagnosis.loc[diagnosis.groupby('RID').AGE.idxmax()].rename(columns={'AGE':'Age_DX_Last', 'DX':'DX_Last'})

# Merge first and last diagnosis
df_visits = df_visits.merge(first_diagnosis[['RID', 'Age_DX_BL','DX_BL']], on='RID', how = 'left')
df_visits = df_visits.merge(last_diagnosis[['RID', 'Age_DX_Last','DX_Last']], on='RID',how = 'left')
replace_first = (df_visits.AGE < df_visits.Age_DX_BL) & (df_visits.DX_BL == 'CN')
df_visits['DX_Extrapolated'] = df_visits['DX']
df_visits.loc[replace_first,'DX_Extrapolated'] = df_visits.loc[replace_first,'DX_BL']
replace_last = (df_visits.AGE > df_visits.Age_DX_Last) & (df_visits.DX_Last == 'AD')
df_visits.loc[replace_last,'DX_Extrapolated'] = df_visits.loc[replace_last,'DX_Last']

# Interpolate Diagnosis
matching_back = 'backward'
Diagnosis_backward = pd.merge_asof(
    left=df_visits.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
    right=df_visits.sort_values('AGE')[['RID','AGE','DX_Extrapolated']].dropna(axis=0, subset=['AGE','DX_Extrapolated']), # variable to find closest match of
    by=['RID'],
    allow_exact_matches=True,
    on='AGE', direction=matching_back,
    ).rename(columns={'DX_Extrapolated' : 'DX_' + matching_back}).sort_values('RID')

matching_forward = 'forward'
Diagnosis_forward = pd.merge_asof(
    left=df_visits.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
    right=df_visits.sort_values('AGE')[['RID','AGE','DX_Extrapolated']].dropna(axis=0, subset=['AGE','DX_Extrapolated']), # variable to find closest match of
    by=['RID'],
    allow_exact_matches=True,
    on='AGE', direction=matching_forward,
        ).rename(columns={'DX_Extrapolated' : 'DX_' + matching_forward }).sort_values('RID')

Diagnosis_backward = Diagnosis_backward.drop_duplicates(subset = ['RID','AGE'],keep = 'first')
Diagnosis_forward = Diagnosis_forward.drop_duplicates(subset = ['RID','AGE'],keep = 'first')

# Merge first and last diagnosis
df_visits = df_visits.merge(Diagnosis_forward.drop(columns=['index']), on=['RID', 'AGE'],how = 'left',validate='m:1')
df_visits = df_visits.merge(Diagnosis_backward.drop(columns=['index']), on=['RID', 'AGE'],how = 'left',validate='m:1')
idx = df_visits['DX_forward']==df_visits['DX_backward']
df_visits['DX_ExIn'] = df_visits['DX_Extrapolated']
df_visits.loc[idx, 'DX_ExIn'] = df_visits.loc[idx, 'DX_forward']

matching_direction = 'nearest'
matching_tolerance = 1.0

# Match-and-merge `Age`
df_visits['AGE_with_DX'] = df_visits['AGE']
df_visits.loc[df_visits['DX_ExIn'].isna(),'AGE_with_DX'] = np.NaN
iD = pd.merge_asof(
        left=df_visits.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
        right=df_visits.sort_values('AGE')[['RID','AGE','DX_ExIn','AGE_with_DX']].dropna(axis=0, subset=['AGE','DX_ExIn']), # variable to find closest match of
        by=['RID'],
        allow_exact_matches=True,
        on='AGE', direction=matching_direction,
        tolerance=matching_tolerance,
                ).rename(columns={'DX_ExIn' : 'DX_' + matching_direction + '_' + str(matching_tolerance)}).sort_values('RID')

iD = iD.drop(columns=['index']).drop_duplicates()
iD.rename(columns={'AGE_with_DX' : 'Years_to_DX'}, inplace=True)
iD['Years_to_DX'] = iD['Years_to_DX'] - iD['AGE']
df_visits = df_visits.merge(iD, on=['RID','AGE'], how='left')

matching_tolerance = 2.0

# Match and Merge FBB
df_visits['AGE_with_FBB'] = df_visits['AGE']
df_visits.loc[df_visits['FBB'].isna(),'AGE_with_FBB'] = np.NaN
iD = pd.merge_asof(
        left=df_visits.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
        right=df_visits.sort_values('AGE')[['RID','AGE','FBB','AGE_with_FBB']].dropna(axis=0, subset=['AGE','FBB']), # variable to find closest match of
        by=['RID'],
        allow_exact_matches=True,
        on='AGE', direction=matching_direction,
        tolerance=matching_tolerance,
                ).rename(columns={'FBB' : 'FBB_' + matching_direction}).sort_values('RID')
iD = iD.drop(columns=['index']).drop_duplicates()
iD.rename(columns={'AGE_with_FBB' : 'Years_to_FBB'}, inplace=True)
iD['Years_to_FBB'] = iD['Years_to_FBB'] - iD['AGE']
df_visits = df_visits.merge(iD, on=['RID','AGE'], how='left')

# Match and Merge AV45
df_visits['AGE_with_AV45'] = df_visits['AGE']
df_visits.loc[df_visits['AV45'].isna(),'AGE_with_AV45'] = np.NaN
iD = pd.merge_asof(
        left=df_visits.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
        right=df_visits.sort_values('AGE')[['RID','AGE','AV45','AGE_with_AV45']].dropna(axis=0, subset=['AGE','AV45']), # variable to find closest match of
        by=['RID'],
        allow_exact_matches=True,
        on='AGE', direction=matching_direction,
        tolerance=matching_tolerance,
                ).rename(columns={'AV45' : 'AV45_' + matching_direction}).sort_values('RID')
iD = iD.drop(columns=['index']).drop_duplicates()
iD.rename(columns={'AGE_with_AV45' : 'Years_to_AV45'}, inplace=True)
iD['Years_to_AV45'] = iD['Years_to_AV45'] - iD['AGE']
df_visits = df_visits.merge(iD, on=['RID','AGE'], how='left')

# Match and Merge PIB
df_visits['AGE_with_PIB'] = df_visits['AGE']
df_visits.loc[df_visits['PIB'].isna(),'AGE_with_PIB'] = np.NaN
iD = pd.merge_asof(
        left=df_visits.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
        right=df_visits.sort_values('AGE')[['RID','AGE','PIB','AGE_with_PIB']].dropna(axis=0, subset=['AGE','PIB']), # variable to find closest match of
        by=['RID'],
        allow_exact_matches=True,
        on='AGE', direction=matching_direction,
        tolerance=matching_tolerance,
                ).rename(columns={'PIB' : 'PIB_' + matching_direction}).sort_values('RID')
iD = iD.drop(columns=['index']).drop_duplicates()
iD.rename(columns={'AGE_with_PIB' : 'Years_to_PIB'}, inplace=True)
iD['Years_to_PIB'] = iD['Years_to_PIB'] - iD['AGE']
df_visits = df_visits.merge(iD, on=['RID','AGE'], how='left')

# Match and Merge ABETA CSF
df_visits['ABETA'] = df_visits['ABETA'].replace('>1700', "1700")
df_visits['ABETA'] = df_visits['ABETA'].str.extract('(\d+)', expand=False).astype(float)
df_visits['AGE_with_ABETA'] = df_visits['AGE']
df_visits.loc[df_visits['ABETA'].isna(),'AGE_with_ABETA'] = np.NaN
iD = pd.merge_asof(
        left=df_visits.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
        right=df_visits.sort_values('AGE')[['RID','AGE','ABETA','AGE_with_ABETA']].dropna(axis=0, subset=['AGE','ABETA']), # variable to find closest match of
        by=['RID'],
        allow_exact_matches=True,
        on='AGE', direction=matching_direction,
        tolerance=matching_tolerance,
                ).rename(columns={'ABETA' : 'ABETA_' + matching_direction}).sort_values('RID')
iD = iD.drop(columns=['index']).drop_duplicates()
iD.rename(columns={'AGE_with_ABETA' : 'Years_to_ABETA'}, inplace=True)
iD['Years_to_ABETA'] = iD['Years_to_ABETA'] - iD['AGE']
df_visits = df_visits.merge(iD, on=['RID','AGE'], how='left')

# Create singular summary column of Amyloid status 
conditions = [
    ((df_visits['FBB_nearest'] <= 1.08) | (df_visits['AV45_nearest'] <= 1.11) | (df_visits['PIB_nearest'] <= 1.5) | (df_visits['ABETA_nearest'] >= 980)),
    ((df_visits['FBB_nearest'] <= 1.08) | (df_visits['AV45_nearest'] >= 1.11) | (df_visits['PIB_nearest'] >= 1.5) | (df_visits['ABETA_nearest'] <= 980)),
]
choices = [0, 1]
df_visits['AMYLOID_SUMMARY'] = np.select(conditions, choices, default=np.NaN)

""" MERGE ALL """
# Merge methylation and cog into visits
data = df_visits.rename(columns={'VISCODE':'VISCODE2'}).merge(dnamCog.drop(columns=['VISCODE','EXAMDATE','PTDOB']),how='outer',on=['RID','VISCODE2'])

# Merge data with imaging
data = data.drop(columns=['FLDSTRENG']).merge(adniLong.rename(columns={'VISCODE':'VISCODE2'}), how='outer', on=['RID','VISCODE2'])

# ffill and bfill to eliminate NAs
data['EXAMDATE'] = data['EXAMDATE'].fillna(data['SCANDATE'])
data['PTDOB'] = data.groupby('RID')['PTDOB'].transform(lambda x: x.ffill().bfill())
data['AGE'] = data['AGE'].fillna((data['EXAMDATE'] - data['PTDOB']).dt.days / 365.25)
data['PTGENDER'] = data.groupby('RID')['PTGENDER'].transform(lambda x: x.ffill().bfill())
data['PTEDUCAT'] = data.groupby('RID')['PTEDUCAT'].transform(lambda x: x.ffill().bfill())
data['APOE4'] = data.groupby('RID')['APOE4'].transform(lambda x: x.ffill().bfill())

# Match and Merge AGE_SAMPLE
data['AGE_with_DNAm'] = data['AGE']
data.loc[data['AGE_SAMPLE'].isna(),'AGE_with_DNAm'] = np.NaN
iD = pd.merge_asof(
        left=data.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
        right=data.sort_values('AGE')[['RID','AGE','AGE_SAMPLE','AGE_with_DNAm']].dropna(axis=0, subset=['AGE','AGE_SAMPLE']), # variable to find closest match of
        by=['RID'],
        allow_exact_matches=True,
        on='AGE', direction=matching_direction, # no matching_tolerance selected
                ).rename(columns={'AGE_SAMPLE' : 'AGE_SAMPLE_' + matching_direction}).sort_values('RID')
iD = iD.drop(columns=['index']).drop_duplicates()
iD.rename(columns={'AGE_with_DNAm' : 'Years_to_DNAm'}, inplace=True)
iD['Years_to_DNAm'] = iD['Years_to_DNAm'] - iD['AGE']
data = data.merge(iD, on=['RID','AGE'], how='left')

# Match and Merge clock ages and age accelerations
clockCols = ['HorvathS2013','HorvathS2013_res',
             'HorvathS2018','HorvathS2018_res',
             'HannumG2013','HannumG2013_res',
             'ShirebyG2020','ShirebyG2020_res',
             'LevineM2018','LevineM2018_res',
             'DunedinPACE',
             'LuA2019','LuA2019_res']

for col in clockCols:
        iD = pd.merge_asof(
                left=data.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
                right=data.sort_values('AGE')[['RID','AGE',col,'AGE_with_DNAm']].dropna(axis=0, subset=['AGE',col]), # variable to find closest match of
                by=['RID'],
                allow_exact_matches=True,
                on='AGE', direction=matching_direction, # no matching_tolerance selected
                        ).rename(columns={col : col + "_" + matching_direction}).sort_values('RID')
        iD = iD.drop(columns=['index']).drop_duplicates()
        data = data.merge(iD.drop(columns=['AGE_with_DNAm']), on=['RID','AGE'], how='left')

data['AVRECTOT'] = data['AVDELTOT'] + 15 - data['AVDELERR2']

# Match and Merge AVLT
matching_tolerance = 1
data['AGE_with_AVLT'] = data['AGE']
data.loc[data['AVRECTOT'].isna(),'AGE_with_AVLT'] = np.NaN
iD = pd.merge_asof(
        left=data.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
        right=data.sort_values('AGE')[['RID','AGE','AVRECTOT','AGE_with_AVLT']].dropna(axis=0, subset=['AGE','AVRECTOT']), # variable to find closest match of
        by=['RID'],
        allow_exact_matches=True,
        on='AGE', direction=matching_direction,
        tolerance=matching_tolerance).rename(columns={'AVRECTOT' : 'AVRECTOT' + "_" + matching_direction}).sort_values('RID')
iD = iD.drop(columns=['index']).drop_duplicates()
iD.rename(columns={'AGE_with_AVLT' : 'Years_to_AVLT'}, inplace=True)
iD['Years_to_AVLT'] = iD['Years_to_AVLT'] - iD['AGE']
data = data.merge(iD, on=['RID','AGE'], how='left')

avlt = ['AVDELTOT','AVDELERR2']
for col in avlt:
        iD = pd.merge_asof(
                left=data.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
                right=data.sort_values('AGE')[['RID','AGE',col,'AGE_with_AVLT']].dropna(axis=0, subset=['AGE',col]), # variable to find closest match of
                by=['RID'],
                allow_exact_matches=True,
                on='AGE', direction=matching_direction, tolerance=matching_tolerance,
                        ).rename(columns={col : col + "_" + matching_direction}).sort_values('RID')
        iD = iD.drop(columns=['index']).drop_duplicates()
        data = data.merge(iD.drop(columns=['AGE_with_AVLT']), on=['RID','AGE'], how='left')

# Match and Merge PHC Cog
matching_tolerance = 1
data['AGE_with_Cog'] = data['AGE']
data.loc[data['PHC_MEM'].isna(),'AGE_with_Cog'] = np.NaN
iD = pd.merge_asof(
        left=data.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
        right=data.sort_values('AGE')[['RID','AGE','PHC_MEM','AGE_with_Cog']].dropna(axis=0, subset=['AGE','PHC_MEM']), # variable to find closest match of
        by=['RID'],
        allow_exact_matches=True,
        on='AGE', direction=matching_direction,
        tolerance=matching_tolerance).rename(columns={'PHC_MEM' : 'PHC_MEM' + "_" + matching_direction}).sort_values('RID')
iD = iD.drop(columns=['index']).drop_duplicates()
iD.rename(columns={'AGE_with_Cog' : 'Years_to_Cog'}, inplace=True)
iD['Years_to_Cog'] = iD['Years_to_Cog'] - iD['AGE']
data = data.merge(iD, on=['RID','AGE'], how='left')

cog = ['PHC_EXF','PHC_LAN', 'PHC_VSP']
for col in cog:
        iD = pd.merge_asof(
                left=data.sort_values('AGE')[['RID','AGE']].dropna(subset=['AGE']).reset_index(),                            # all entries
                right=data.sort_values('AGE')[['RID','AGE',col,'AGE_with_Cog']].dropna(axis=0, subset=['AGE',col]), # variable to find closest match of
                by=['RID'],
                allow_exact_matches=True,
                on='AGE', direction=matching_direction, tolerance=matching_tolerance,
                        ).rename(columns={col : col + "_" + matching_direction}).sort_values('RID')
        iD = iD.drop(columns=['index']).drop_duplicates()
        data = data.merge(iD.drop(columns=['AGE_with_Cog']), on=['RID','AGE'], how='left') 

""" SORT """
data = data.sort_values(['RID','AGE'], ascending=[True, True]).reset_index(drop=True)

""" SAVE """
# Save CSV
data.to_csv(dir + 'data/ADNI/MRI_DNA_COG_LONG.csv')
