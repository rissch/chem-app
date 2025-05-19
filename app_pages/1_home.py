import streamlit as st
import pandas as pd


# set page config
st.set_page_config(
    page_title="Chemistry Hub",
    page_icon=":test_tube:",
    layout="wide",
    initial_sidebar_state="expanded"
)

def display_home():
    st.title("Welcome to Chemistry Hub")
    st.write("Explore various tools to simplify complex chemical calculations.")
    
    st.markdown("[Educational Tools](edu_tools)")
    st.markdown("- Periodic Table")
    st.markdown("- Practice Chemical Equations")
    st.markdown("- Visual Molecule Builder")
    st.markdown("- Flashcards")


    st.markdown("[Learning Hub](learning_hub)")
    st.markdown("- Chemistry Concepts")
    st.markdown("- Problem Solving Guides")
    st.markdown("- Downloadable Reference Materials")
    st.markdown("- Visuals and tutorials")
    

    st.markdown("[Basic Calculators](basic_calculators)")
    st.markdown("- Molar Mass Calculator")
    st.markdown("- ⚖️ Stoichiometry Calculator")
    st.markdown("- 💧 Molarity Calculator")
    st.markdown("- 🌡️ Molality Calculator")
    st.markdown("- 📏 Normality Calculator")
    st.markdown("- 🌊 PPM Calculator")


    st.markdown("[Analytical Tools](analytical_tools)")
    st.markdown("- 📡 Spectroscopy Tool")
    st.markdown("- 🧪 Chromatography Calculator")
    st.markdown("- ⚗️ Titration Calculator")
    st.markdown("- ⚗️ pH Calculator")


    st.markdown("[Physical Chemistry](physical_chemistry)")
    st.markdown("- 🌬️ Ideal Gas Law")
    st.markdown("- 🔥 Thermodynamics")
    st.markdown("- ⏱️ Kinetics")
    st.markdown("- ⚖️ Equilibrium")
    st.markdown("- 🔋 Electrochemistry")
    st.markdown("- 🔬 Quantum Mechanics")


    st.markdown("[Reaction Tools](reaction_tools)")
    st.markdown("- 🔀 Equation Balancer")
    st.markdown("- ⚖️ Limiting Reagent")
    st.markdown("- 🔥 Reaction Enthalpy")
    st.markdown("- 📈 Reaction Rate")
    st.markdown("- 🎛️ Equilibrium Constant")

    st.markdown("[Unit Converters](unit_converters)")
    st.markdown("- 🧪 Mass")
    st.markdown("- 🧴 Volume")
    st.markdown("- 🌡️ Temperature")
    st.markdown("- 🌬️ Pressure")
    st.markdown("⚡ Energy")

    st.markdown("[Blogs](blogs)")
    st.markdown("- 📚 Study Tips")
    st.markdown("- 📝 Exam Preparation")
    st.markdown("- 🧯 Lab Safety")
    st.markdown("- 💼 Career Guides")
    
    st.markdown("[About us](about_us)")
    st.markdown("[Contact us](contact_us)")

    st.markdown("###")

    # Footer HTML & CSS
    footer = """
    <style>
    /* Remove default padding/margin */
    .reportview-container .main {
        padding-bottom: 0px;
    }

    footer {
        visibility: hidden;
    }

    .footer {
        position: fixed;
        bottom: 0;
        left: 0;
        right: 0;
        background-color: #2d2d2d;  /* neutral dark */
        color: #ddd;  /* light gray text */
        text-align: center;
        padding: 12px 0;
        font-size: 14px;
        z-index: 1000;
        border-top: 1px solid #444;
    }
    .footer a {
        color: #aaa;
        text-decoration: none;
        margin: 0 12px;
    }
    .footer a:hover {
        color: white;
        text-decoration: underline;
    }
    </style>

    <div class="footer">
        🔬 Built by <strong>Molaritycalculator.net</strong> |
        <a href="/contact_us">Contact</a> |
        <a href="[https://github.com/your-repo]" target="_blank">GitHub</a>
    </div>
    """

    st.markdown(footer, unsafe_allow_html=True)


if __name__ == "__main__":
    display_home()
