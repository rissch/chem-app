import json
import streamlit as st
from streamlit_option_menu import option_menu
import chempy
from chempy import balance_stoichiometry
from chempy.util import periodic

def display_reaction_tools():
    calc = option_menu(
        menu_title="Reaction Tools",
        options = [
            "🔀 Equation Balancer",
            "⚖️ Limiting Reagent",
            "🔥 Reaction Enthalpy",
            "📈 Reaction Rate",
            "🎛️ Equilibrium Constant"
        ],
        icons = [
            "none",      
            "none",
            "none",        
            "none",  
            "none"           
        ],
        orientation="horizontal"
    )

    if calc == "🔀 Equation Balancer":
    

        st.subheader("⚖️ Chemical Equation Balancer")

        st.markdown("Enter your **unbalanced chemical equation** below. Use `+` to separate reactants/products, and `->` to separate sides.")
        
        with open("ChemicalEquations.json", "r") as f:
            equations = json.load(f)[0]['balance_reaction_examples']
        
        option = st.checkbox("Enter custom equation", value=False)
        if option:
            user_eq = st.text_input("Enter your equation", key="custom_eq")
            st.code(f"Example : H2 + O2 -> H2O", language="markdown")
        else:
            user_eq = st.selectbox("Select Example Equation", equations.keys())
            st.code(f"Selected : {user_eq}", language="markdown")


        if user_eq:
            try:
                # Parse equation
                if '->' not in user_eq:
                    st.error("❌ Equation must contain '->' to separate reactants and products.")
                else:
                    reactants_str, products_str = user_eq.split('->')
                    reactants = [r.strip() for r in reactants_str.split('+')]
                    products = [p.strip() for p in products_str.split('+')]

                    reac_set, prod_set = set(reactants), set(products)

                    # Balance using chempy
                    reac_coeffs, prod_coeffs = balance_stoichiometry(reac_set, prod_set)

                    # Display result
                    def format_side(side, coeffs):
                        return ' + '.join([f"{coeffs[chem]} {chem}" if coeffs[chem] != 1 else chem for chem in side])

                    balanced_eq = f"{format_side(reactants, reac_coeffs)} -> {format_side(products, prod_coeffs)}"
                    st.success("✅ Balanced Equation:")
                    st.code(balanced_eq)

            except Exception as e:
                st.error(f"⚠️ Unable to balance the equation. Reason: {e}")

    elif calc == "⚖️ Limiting Reagent":
        st.subheader("⚖️ Limiting Reagent Calculator")
        reactant1_moles = st.number_input("Moles of Reactant 1", min_value=0.0)
        reactant1_coeff = st.number_input("Coefficient of Reactant 1", min_value=1.0)
        reactant2_moles = st.number_input("Moles of Reactant 2", min_value=0.0)
        reactant2_coeff = st.number_input("Coefficient of Reactant 2", min_value=1.0)

        ratio1 = reactant1_moles / reactant1_coeff
        ratio2 = reactant2_moles / reactant2_coeff

        if ratio1 < ratio2:
            st.success("Reactant 1 is the limiting reagent.")
        elif ratio2 < ratio1:
            st.success("Reactant 2 is the limiting reagent.")
        else:
            st.success("Both reactants are present in stoichiometric amounts.")

    elif calc == "🔥 Reaction Enthalpy":
        st.subheader("🔥 Reaction Enthalpy Calculator")
        st.latex("\Delta H = \Sigma H_{products} - \Sigma H_{reactants}")
        product_enthalpy = st.number_input("Total enthalpy of products (kJ/mol)")
        reactant_enthalpy = st.number_input("Total enthalpy of reactants (kJ/mol)")
        delta_H = product_enthalpy - reactant_enthalpy
        st.success(f"ΔH = {delta_H:.3f} kJ/mol")

    elif calc == "📈 Reaction Rate":
        st.subheader("⚗️ Reaction Rate Calculator")
        st.latex("rate = \Delta [A]/\Delta t")
        delta_conc = st.number_input("Change in concentration Δ[A] (mol/L)")
        delta_time = st.number_input("Change in time Δt (s)")

        if delta_time == 0:
            st.error("Time interval cannot be zero.")
        else:
            rate = delta_conc / delta_time
            st.success(f"Rate = {rate:.3f} mol/L·s")

    elif calc == "🎛️ Equilibrium Constant":
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