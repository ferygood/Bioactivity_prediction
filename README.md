# Bioactivity_prediction

This project is to build a machine learning model using the [ChEMBL](https://www.ebi.ac.uk/chembl/) bioactivity database. ChEMBL is a manually curated database of bioactive molecules with drug-like properties. It brings together chemical, bioactivity and genomic data to aid the translation of genomic information into effective new drugs. Current version is ChEMBL 32.

## Download Bioactivity Data
Use their python package [chembl_webresource_client](https://github.com/chembl/chembl_webresource_client), in the first example, we will search target for Acetycholinesterase. The script for first part is in `download_data.py`. It downloads the data and first selects the information including **molecule_chembl_id**, **canonical_smiles**, **standard_value**. Then based on the value of standard value, classifying the bioactivity class into inactive(>=1000), active(<=1000), and intermediate (others). Finally, save the dataframe into `data/chembl220_bioactivity.csv`. 

## Exploratory Data Analysis (Chemical Space Analysis)
Apply Lipinski rule of 5 on data. Normalize IC50 to pIC50. 

## Containerization
Use docker to contain ... this app