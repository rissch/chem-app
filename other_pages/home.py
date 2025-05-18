import streamlit as st
import pandas as pd

from utils import filter_compounds_with_type, load_compounds


def display_home():
    st.title("Welcome to Chemistry Hub")
    st.write("Explore various tools to simplify complex chemical calculations.")

   
