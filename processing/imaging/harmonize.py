import pandas as pd
import numpy as np
from scipy import stats
from neuroHarmonize import harmonizationLearn, harmonizationApply

# Read data
dir = '~/Projects/CoMeT/'
adni = pd.read_csv(dir + 'data/ADNI/MRI_DNA_COG_XS.csv',dtype={'RID':int, 'VISCODE':str}, parse_dates=['EXAMDATE'],low_memory=False)

# Encode categorical variables 
adni['SITE'] = adni['FLDSTRENG'].astype('str')
adni['PTGENDER_M'] = pd.factorize(adni['PTGENDER'])[0]
adni['DX_nearest_1.0_ordinal'] = pd.factorize(adni['DX_nearest_1.0'])[0] + 1

# Extract FreeSurfer ROIs that match the regex pattern '^ST\\d+TA$'
names = adni.columns.tolist()
rois = [name for name in names if pd.Series(name).str.match(r'^ST\d+(?:TA)$').any()]
rois = rois + ['ST29SV','ST88SV','ST12SV','ST71SV','ST10CV']

# ROIs to exclude due to large amounts of missing data
rois_missing = ["ST123TA", "ST22TA", "ST64TA", "ST81TA"]

# Remove missing ROIs from the list
rois = [roi for roi in rois if roi not in rois_missing]
rois_combat = [roi + '_H' for roi in rois]

# Data to harmonize
cn = adni.loc[adni['Case']==0].reset_index(drop=True)
case = adni.loc[adni['Case']==1].reset_index(drop=True)

# Controls
data = cn[rois]
data = data.to_numpy()

# Specifying the batch variable as well as biological covariate to preserve:
covars = cn[['SITE','AGE','PTGENDER_M']] 
covars = pd.DataFrame(covars)  

# Run harmonization
model, cn_h = harmonizationLearn(data, covars)
cn_h = pd.DataFrame(cn_h,columns=rois_combat)
cn_h['RID'] = cn['RID']

# Prepare held out data 
data = case[rois].to_numpy()
covars = case[['SITE','AGE','PTGENDER_M']]

# Apply harmonization
case_h = harmonizationApply(data, covars, model)
case_h = pd.DataFrame(case_h,columns=rois_combat)
case_h['RID'] = case['RID']

# Merge
cn_combat = pd.merge(cn,cn_h,on='RID',how='inner')
case_combat = pd.merge(case,case_h,on='RID',how='inner')
adni_combat = pd.concat([cn_combat,case_combat],axis=0)

# Visualize harmonization
import matplotlib.pyplot as plt
import seaborn as sns
plt.figure()
site_01 = stats.norm.rvs(size=10000, loc=model['gamma_bar'][0], scale=np.sqrt(model['t2'][0]))
site_02 = stats.norm.rvs(size=10000, loc=model['gamma_bar'][1], scale=np.sqrt(model['t2'][1]))
sns.kdeplot(site_01,color='orange', label='1.5T Prior')
sns.kdeplot(model['gamma_hat'][0, :], color='orange', label='1.5T Observed', linestyle='--')
sns.kdeplot(site_02,color='brown', label='3T Prior')
sns.kdeplot(model['gamma_hat'][1, :], color='brown', label='3T Observed', linestyle='--')
plt.legend()
plt.title('Prior and Observed Batch Effects in ADNI MRI')

# Save to PDF
plt.savefig("harmonization.pdf", format="pdf")
plt.show()

# (Optional) Close the figure
plt.close()

# Write
adni_combat.to_csv(dir + 'data/ADNI/MRI_DNA_COG_XSH.csv')