import pandas as pd
from scipy.stats import pearsonr

print("reading..")
fname ="chr7training.RandOncontacts.gz.False.11.1500000.10001.1.250000.cont_with_CTCF341430.5000"
fname = "chr7training.RandOncontacts.gz.False.11.1500000.10001.1.25000.cont_with_CTCF341430.5000"
contact_data = pd.read_csv("/mnt/scratch/ws/psbelokopytova/202001051010polina_data/3DPredictor/out/K562/5KB/all_predictors/chr7training.RandOncontacts.gz.False.11.1500000.10001.1.250000.cont_with_CTCF341430.5000", sep="\t")
print(contact_data.keys())
dists = [25000, 75000, 100000]
corrs = []
n_contacts_list=[]
for dist in dists:
    print(dist)
    data_query = contact_data[contact_data["contact_dist"]==dist]
    # data_query=contact_data
    correlation = pearsonr(data_query["contact_count"], data_query["CTCF_W"])
    corrs.append(correlation)
    n_contacts_list.append(len(data_query))
corr_df = pd.DataFrame({'chr':["chr14"]*len(dists), 'dist':dists, 'n_contacts':n_contacts_list, 'pearsonR':corrs})
corr_df.to_csv("correlations_K562.txt", index=False, sep="\t")
# print(correlation)
