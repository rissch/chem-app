# imports
import streamlit as st
from streamlit_option_menu import option_menu

from other_pages.learning_hub import display_learning_hub

# Set page config
st.set_page_config(page_title="Chemistry Hub", layout="wide")

from other_pages.analytics_tools import display_analytics_tools
from other_pages.basic_calculator import display_basic_calculator
from other_pages.edu_tools import display_edu_tools
from other_pages.home import display_home
from other_pages.physical_chemistry import display_physical_chemistry
from other_pages.reaction_tools import display_reaction_tools
from other_pages.unit_converters import display_unit_converters
from other_pages.blog import display_blogs

# Sidebar Navigation
with st.sidebar:
    selected = option_menu(
        menu_title="Chemistry Hub",
        options=[
            "Home",
            "Educational Tools",
            "Learning Hub",
            "Basic Calculators",
            "Analytical Tools",
            "Physical Chemistry",
            "Reaction Tools",
            "Unit Converters",
            "Blogs",
            "About Us",
            "Contact",            
        ],
        icons=[
        "house", "book", "journal-text", "calculator", "graph-up", "beaker", "activity",
        "sliders", "file-earmark-text", "info-circle","envelope",    
    ],
        menu_icon="chemistry",
        default_index=0
    )

# Home Page
if selected == "Home":
    display_home()

elif selected == "Educational Tools":
    display_edu_tools()

elif selected == "Learning Hub":
    display_learning_hub()

# Basic Calculators
elif selected == "Basic Calculators":
    display_basic_calculator()

# Analytical Tools
elif selected == "Analytical Tools":
    display_analytics_tools()

# Physical Chemistry
elif selected == "Physical Chemistry":
    display_physical_chemistry()

# Reaction Tools
elif selected == "Reaction Tools":
    display_reaction_tools()

# Unit Converters
elif selected == "Unit Converters":
    display_unit_converters()

# blogs
elif selected == "Blogs":
    display_blogs()

# About Us
elif selected == "About Us":
    st.title("About Chemistry Calculator Hub")
    st.write("This platform is designed to help learners, educators, and researchers...")

# Contact
elif selected == "Contact":
    st.title("Contact Us")
    st.write("For feedback, collaborations, or support, email us at: support@chemcalc.com")


