import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt

df = pd.read_csv("data/chembl220_bioactivity.csv")

# create lipinski function
# Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation
# https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five
def lipinski(smiles, verbose=False):
    
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1,1)
    i = 0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
        
        row = np.array([
            desc_MolWt,
            desc_MolLogP,
            desc_NumHDonors,
            desc_NumHAcceptors
        ])

        if(i==0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i+=1

    columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors 
      
df_lipinski = lipinski(df.canonical_smiles)

# combine original dataframe with lipinski value
df_combine = pd.concat([df, df_lipinski], axis=1)
print(df_combine)


# To make distribution more uniformly distributed, convert IC50 to pIC50 (-log10(IC50))
# inspired: https://github.com/chaninlab/estrogen-receptor-alpha-qsar/blob/master/02_ER_alpha_RO5.ipynb
def pIC50(input):
    pIC50 = []
    for i in input['standard_value_norm']:
        molar = i*(10**-9) # converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
    
    return x

def norm_value(input):
    norm = []
    for i in input['standard_value']:
        if i > 100000000:
            i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)

    return x 

# first normalized the data, standard_value greater than 100000000 will be fixed at 10**8 otherwise 
# the negative logarithmic value will become negative
df_norm = norm_value(df_combine)
df_norm_pIC50 = pIC50(df_norm)

# removing intermediate class
df_2class = df_norm_pIC50[df_norm_pIC50['class'] != 'intermediate']

# visualize boxplot: frequency plot
plt.figure(figsize=(5.5, 5.5))
sns.countplot(x='class', data=df_2class, edgecolor='black')
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')
plt.savefig('images/plot_bioactivity_class.pdf')

# boxenplot pIC50
plt.figure(figsize=(5.5, 5.5))
sns.boxenplot(x = 'class', y = 'pIC50', data = df_2class)
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')
plt.savefig('images/plot_ic50.pdf')

# Statistical analysis: Mann-Whitney U Test
def mannwhitney(descriptor, verbose=False):
  # https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python/
  from numpy.random import seed
  from numpy.random import randn
  from scipy.stats import mannwhitneyu

# seed the random number generator
  seed(1)

# actives and inactives
  selection = [descriptor, 'class']
  df = df_2class[selection]
  active = df[df['class'] == 'active']
  active = active[descriptor]

  selection = [descriptor, 'class']
  df = df_2class[selection]
  inactive = df[df['class'] == 'inactive']
  inactive = inactive[descriptor]

# compare samples
  stat, p = mannwhitneyu(active, inactive)
  #print('Statistics=%.3f, p=%.3f' % (stat, p))

# interpret
  alpha = 0.05
  if p > alpha:
    interpretation = 'Same distribution (fail to reject H0)'
  else:
    interpretation = 'Different distribution (reject H0)'
  
  results = pd.DataFrame({'Descriptor':descriptor,
                          'Statistics':stat,
                          'p':p,
                          'alpha':alpha,
                          'Interpretation':interpretation}, index=[0])
  filename = 'mannwhitneyu_' + descriptor + '.csv'
  results.to_csv(filename)

  return results

print(mannwhitney('MW'))

# Check if all of the 4 Lipinski's descriptors exhibited statistically difference between the actives and inactives.

df_2class.to_csv('data/chembl220_bioactivity_2class_pIC50.csv')