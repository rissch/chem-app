import pybel
import sqlite3
import pandas as pd
from rdkit import Chem
import streamlit as st
from rdkit.Chem import GetPeriodicTable


@st.cache_data
def load_compounds():
    conn = sqlite3.connect("chemistry.db")
    df = pd.read_sql_query("SELECT Formula, SMILES FROM compounds", conn)
    conn.close()
    df["mol"] = df["SMILES"].apply(Chem.MolFromSmiles)
    return df

@st.cache_data
def load_formulas_only():
    conn = sqlite3.connect("chemistry.db")
    df = pd.read_sql_query("SELECT Formula FROM compounds", conn)
    conn.close()
    return df["Formula"].tolist()

@st.cache_data
def get_smiles(formula):
    conn = sqlite3.connect("chemistry.db")
    query = f"SELECT SMILES FROM compounds WHERE Formula = ?"
    result = pd.read_sql_query(query, conn, params=(formula,))
    conn.close()
    return result["SMILES"].iloc[0] if not result.empty else None

@st.cache_data
def filter_compounds_with_type(type):
    conn = sqlite3.connect("chemistry.db")
    df = pd.read_sql_query(f"SELECT Formula, SMILES FROM compounds WHERE Type = '{type}'", conn)
    conn.close()
    df["mol"] = df["SMILES"].apply(Chem.MolFromSmiles)
    return df

@st.cache_data
def get_atoms():
    ptable = GetPeriodicTable()

    # Get all elements up to atomic number 118
    atom_list = []
    for atomic_num in range(1, 119):  # 1 to 118 (Hydrogen to Oganesson)
        symbol = ptable.GetElementSymbol(atomic_num)
        name = ptable.GetElementName(atomic_num)
        atom_list.append((symbol, name))

    # Print the list
    for symbol, name in atom_list:
        print(f"{symbol}: {name}")

    return atom_list

