import streamlit as st
import pandas as pd
import numpy as np

df = pd.read_csv("data/chembl220_bioactivity.csv")

st.dataframe(df)