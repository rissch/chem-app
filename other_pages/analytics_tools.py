import streamlit as st
from streamlit_option_menu import option_menu
import math

def display_analytics_tools():
    st.title("🔬 Analytical Tools")
    st.write("Explore various analytical tools to simplify complex chemical calculations and data interpretation.")

    calc = option_menu(
        menu_title="Analytical Tools",
        options=["Spectroscopy", "Chromatography", "Titration", "pH"],
        icons=["activity", "bar-chart", "droplet", "thermometer-half"],
        orientation="horizontal"
    )

    if calc == "Spectroscopy":
        st.subheader("📡 Spectroscopy Tool")
        st.write("Estimate absorbance or concentration using **Beer-Lambert Law**:")
        st.latex(r"A = \varepsilon \cdot c \cdot l")
        
        ε = st.number_input("Molar Absorptivity (ε) in L/mol·cm", min_value=0.0, step=0.01)
        c = st.number_input("Concentration (c) in mol/L", min_value=0.0, step=0.01)
        l = st.number_input("Path Length (l) in cm", value=1.0, min_value=0.1, step=0.1)

        if st.button("Calculate Absorbance"):
            A = ε * c * l
            st.success(f"Absorbance (A) = {A:.4f}")

    elif calc == "Chromatography":
        st.subheader("🧪 Chromatography Calculator")
        st.write("Analyze retention time and resolution between compounds.")
        
        rt1 = st.number_input("Retention Time of Compound 1 (min)", min_value=0.0, step=0.1)
        rt2 = st.number_input("Retention Time of Compound 2 (min)", min_value=0.0, step=0.1)
        w1 = st.number_input("Peak Width of Compound 1 (min)", min_value=0.0, step=0.1)
        w2 = st.number_input("Peak Width of Compound 2 (min)", min_value=0.0, step=0.1)

        if st.button("Calculate Resolution"):
            try:
                Rs = (2 * abs(rt2 - rt1)) / (w1 + w2)
                st.success(f"Resolution (Rs) = {Rs:.3f}")
            except ZeroDivisionError:
                st.error("Peak widths cannot be zero.")

    elif calc == "Titration":
        st.subheader("⚗️ Titration Calculator")
        st.write("Calculate unknown concentration from titration data.")

        V1 = st.number_input("Volume of Titrant (mL)", min_value=0.0)
        C1 = st.number_input("Concentration of Titrant (mol/L)", min_value=0.0)
        V2 = st.number_input("Volume of Analyte (mL)", min_value=0.0)

        if st.button("Calculate Analyte Concentration"):
            try:
                C2 = (V1 * C1) / V2
                st.success(f"Concentration of Analyte = {C2:.4f} mol/L")
            except ZeroDivisionError:
                st.error("Analyte volume cannot be zero.")

    elif calc == "pH":
        st.subheader("🧪 pH Calculator")
        st.write("Calculate pH from **[H⁺]** or **[OH⁻]** concentration.")

        mode = st.radio("Choose input type:", ["[H⁺] concentration", "[OH⁻] concentration"])

        if mode == "[H⁺] concentration":
            H = st.number_input("[H⁺] in mol/L", min_value=1e-14, format="%.2e")
            pH = -math.log10(H)  # Correctly using base-10 logarithm
            st.success(f"pH = {pH:.3f}")
        else:
            OH = st.number_input("[OH⁻] in mol/L", min_value=1e-14, format="%.2e")
            pOH = -math.log10(OH)  # Correctly using base-10 logarithm
            pH = 14 - pOH
            st.success(f"pH = {pH:.3f}")
