import streamlit as st
from streamlit_option_menu import option_menu


def display_physical_chemistry():
    calc = option_menu(
        menu_title="Physical Chemistry",
        options=["Ideal Gas Law", "Thermodynamics", "Kinetics", "Equilibrium", "Electrochemistry", "Quantum Mechanics"],
        icons=["wind", "sun"],
        orientation="horizontal"
    )

    st.header(f"{calc}")

    if calc == "Ideal Gas Law":
        st.subheader("🌡️ Ideal Gas Law Calculator")
        st.latex("PV = nRT")
        P = st.number_input("Pressure (P) in atm", min_value=0.0, value=1.0)
        V = st.number_input("Volume (V) in L", min_value=0.0, value=1.0)
        n = st.number_input("Amount of substance (n) in mol", min_value=0.0, value=1.0)
        R = 0.0821  # Ideal gas constant in L·atm/(mol·K)
        T = st.number_input("Temperature (T) in K", min_value=0.0, value=273.15)

        calc_mode = st.radio("Solve for:", ["Pressure", "Volume", "Moles", "Temperature"])

        if calc_mode == "Pressure":
            P = (n * R * T) / V
            st.success(f"Pressure = {P:.3f} atm")
        elif calc_mode == "Volume":
            V = (n * R * T) / P
            st.success(f"Volume = {V:.3f} L")
        elif calc_mode == "Moles":
            n = (P * V) / (R * T)
            st.success(f"Moles = {n:.3f} mol")
        elif calc_mode == "Temperature":
            T = (P * V) / (n * R)
            st.success(f"Temperature = {T:.3f} K")

    elif calc == "Thermodynamics":
        st.subheader("🔥 Thermodynamics Calculator")
        st.latex("\Delta G = \Delta H - T\Delta S")
        H = st.number_input("ΔH (Enthalpy) in kJ/mol", value=0.0)
        S = st.number_input("ΔS (Entropy) in J/mol·K", value=0.0)
        T = st.number_input("Temperature (T) in K", value=298.15)

        G = H - (T * S / 1000)  # Convert S to kJ/mol·K
        st.success(f"ΔG = {G:.3f} kJ/mol")

    elif calc == "Kinetics":
        st.subheader("⚗️ Chemical Kinetics Calculator")
        st.latex("rate = k[A]^n")
        k = st.number_input("Rate constant (k)", value=1.0)
        A = st.number_input("Concentration of A [A] in mol/L", value=1.0)
        n = st.number_input("Order of reaction (n)", value=1.0)

        rate = k * (A ** n)
        st.success(f"Rate = {rate:.3f} mol/L·s")

    elif calc == "Equilibrium":
        st.subheader("⚖️ Equilibrium Constant Calculator")
        st.latex("K = [products]^{coeff} / [reactants]^{coeff}")
        prod_conc = st.number_input("Product concentration (mol/L)", value=1.0)
        prod_coeff = st.number_input("Product coefficient", value=1.0)
        reac_conc = st.number_input("Reactant concentration (mol/L)", value=1.0)
        reac_coeff = st.number_input("Reactant coefficient", value=1.0)

        if reac_conc ** reac_coeff == 0:
            st.error("Reactant concentration raised to its coefficient cannot be zero.")
        else:
            K = (prod_conc ** prod_coeff) / (reac_conc ** reac_coeff)
            st.success(f"Equilibrium constant K = {K:.3f}")

    elif calc == "Electrochemistry":
        st.subheader("🔋 Electrochemistry Calculator")
        st.latex("\Delta G = -nFE")
        n = st.number_input("Number of electrons transferred (n)", value=1.0)
        F = 96485  # Faraday constant in C/mol
        E = st.number_input("Cell potential (E) in volts", value=1.0)

        G = -n * F * E / 1000  # Convert to kJ
        st.success(f"ΔG = {G:.3f} kJ/mol")

    elif calc == "Quantum Mechanics":
        st.subheader("🌀 Quantum Energy Calculator")
        st.latex("E = h\nu")
        h = 6.626e-34  # Planck constant in J·s
        freq = st.number_input("Frequency (ν) in Hz", value=1e14)

        E = h * freq
        st.success(f"Energy = {E:.3e} J")