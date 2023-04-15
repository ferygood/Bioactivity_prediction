import pandas as pd
from chembl_webresource_client.new_client import new_client

# target search for coronavirus
target = new_client.target
target_query = target.search('acetylcholinesterase')
targets = pd.DataFrame.from_dict(target_query)

# select and retrieve bioactivity data for Human Acetylcholinesterase (ID: CHEMBL220)
selected_target = targets.target_chembl_id[0]
#print(selected_target)

# filter activity data
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)

# handling missing data
# if any compounds has missing value for standard_value and canonical_smiles column then drop it
df = df[df.standard_value.notna()]
df = df[df.canonical_smiles.notna()]
df = df.drop_duplicates(['canonical_smiles']) # remove duplicate value
df = df.reset_index(drop=True)

# combine 3 columns including molecule_chembl_id, canonical_smiles, standard_value with bioactivity_class into a dataframe
# first extract these 3 columns
selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df = df[selection]
# labeling bioactivity compound
bioactivity_threshold = []
for i in df.standard_value:
    if float(i) >= 1000:
        bioactivity_threshold.append("inactive")
    elif float(i) <= 1000:
        bioactivity_threshold.append("active")
    else:
        bioactivity_threshold.append("intermediate")

bioactivity_class = pd.Series(bioactivity_threshold, name='class')
df_combine = pd.concat([df, bioactivity_class], axis=1)

# save the bioactivity data to a csv file chembl220_bioactivity.csv
df_combine.to_csv("data/chembl220_bioactivity.csv", index=False)