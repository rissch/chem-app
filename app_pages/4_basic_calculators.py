import json
import streamlit as st
import re

from periodictable import elements
from rdkit.Chem import Draw
from streamlit_option_menu import option_menu
from chempy import Reaction

from utils import filter_compounds_with_type, load_compounds


df = load_compounds()

formulas = df['Formula'].tolist()
# df["mol"] = df["SMILES"].apply(Chem.MolFromSmiles)


def display_basic_calculator():
    calc = option_menu(
        menu_title="Basic Calculators",
        options=[
        "🧮 Molar Mass",
        "⚖️ Stoichiometry",
        "💧 Molarity",
        "🌡️ Molality",
        "📏 Normality",
        "🌊 PPM"
    ],
    icons=["none", "none", "none", "none", "none", "none"],
        orientation="horizontal"
    )
    st.header(f"{calc} Calculator")
    # Call corresponding function later here

    if calc == "🧮 Molar Mass":
        molar_mass_calculator()
    elif calc == "⚖️ Stoichiometry":
        stoichiometry_calculator()
    elif calc == "💧 Molarity":
        molarity_calculator()
    elif calc == "🌡️ Molality":
        molality_calculator()
    elif calc == "📏 Normality":
        normality_calculator()
    elif calc == "🌊 PPM":
        ppm_calculator()

def molar_mass_calculator():
    with st.expander("📘 What is Molar Mass?"):
                st.write("""
                Molar mass is the mass of one mole of a substance. It's calculated by summing the atomic masses
                of all atoms in the formula. This value is used to convert between mass and moles in chemical calculations.
                """)
    st.markdown("Enter a **chemical formula** (e.g., `H2O`, `NaCl`, `C6H12O6`) to calculate its molar mass.")

    formula_type = st.selectbox("Select chemical formula type", ["Organic", "Inorganic", "Acid"], index=0, key="formula_type")

    filtered_compounds = filter_compounds_with_type(formula_type)

    formula = st.selectbox("Select a chemical formula or write yourself", options=filtered_compounds, index=0)
    # st.text_input("Chemical Formula", value="H2O")
    
    def parse_formula(formula):
        # Match elements with optional numbers (e.g., H2, O)
        pattern = r'([A-Z][a-z]?)(\d*)'
        parsed = re.findall(pattern, formula)
        return [(elem, int(count) if count else 1) for elem, count in parsed]

    if formula:
        try:
            parsed = parse_formula(formula)
            total_mass = 0
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("### 🧮 Element Breakdown")
                for symbol, count in parsed:
                    element = getattr(elements, symbol, None)
                    if not element:
                        st.error(f"Element '{symbol}' not found in periodic table.")
                        return
                    atomic_mass = element.mass
                    subtotal = atomic_mass * count
                    total_mass += subtotal
                    st.markdown(f"- **{symbol}**: {count} × {atomic_mass:.3f} g/mol = `{subtotal:.3f}` g/mol")

                st.success(f"**Total Molar Mass of {formula} = {total_mass:.3f} g/mol**")


            # Displaying the structure of the compound
            with col2:
                mol = df[df["Formula"] == formula]["mol"].iloc[0]
                if mol:
                    st.markdown(f"### ⌬ Structure of {formula}")
                    st.image(Draw.MolToImage(mol), caption=f"Structure of {formula}")
                else:
                    st.info("No structure available.")

            

        except Exception as e:
            st.error(f"Error parsing formula: {e}")

def stoichiometry_calculator():
    with st.expander("📘 What is Stoichiometry?"):
        st.write("""
        Stoichiometry is the calculation of reactants and products in chemical reactions. It uses the coefficients
        from a balanced equation to determine the amount of substances involved. This is essential for predicting
        yields and understanding reaction mechanisms.
        """)

    st.write("Enter a balanced chemical equation (e.g., `2 H2 + O2 -> 2 H2O`)")
    # equation_input = st.text_input("Balanced Equation")

    if "selected_formula" not in st.session_state:
        st.session_state.selected_formula = ""

    with open("ChemicalEquations.json", "r") as f:
        eqs_list = json.load(f)[0]['formula_examples']

    # User input (use full string for now)
    # full_input = st.text_input("🔍 Type a compound or formula", value=st.session_state.selected_formula)
    full_input = st.selectbox("🔍 Type a compound or formula", options=eqs_list)

    # # Extract the last word typed
    # last_term = full_input.strip().split()[-1] if full_input.strip() else ""

    # if full_input != st.session_state.selected_formula:
    #     st.session_state.selected_formula = full_input

    # # Only search if there's a last term
    # if last_term:
    #     matches_startswith = df[
    #         df['Formula'].str.startswith(last_term)
    #     ].head(3)

    #     matches_contain = df[
    #         df['Formula'].str.contains(last_term, case=False)
    #     ].head(3)

    #     matches = pd.concat([matches_startswith, matches_contain]).drop_duplicates().reset_index(drop=True)

    #     if not matches.empty:
    #         cols = st.columns(len(matches))  # Create as many columns as there are suggestions

    #         for i, row in matches.iterrows():
    #             label = f"**{row['Formula']}**"
    #             with cols[i]:
    #                 if st.button(label, key=f"suggestion_{i}"):
    #                     tokens = full_input.strip().split()
    #                     tokens[-1] = row['Formula']
    #                     new_input = " ".join(tokens)
    #                     st.session_state.selected_formula = new_input
    #                     st.rerun()

    #     else:
    #         st.info(f"No matches found for '{last_term}'")


    # if st.button("Apply"):
    #     st.session_state.show_stoich = True 

    # if st.session_state.get("show_stoich", False):
    try:
        # Parse the equation
        reaction = Reaction.from_string(full_input)

        # Combine all compounds from reactants and products
        all_species = list(reaction.reac.keys()) + list(reaction.prod.keys())

        if all_species:
            known_substance = st.selectbox("Select known substance", all_species)
            known_moles = st.number_input("Amount in moles of known substance", min_value=0.0, format="%.4f")
            desired_substance = st.selectbox("Select desired substance", all_species)

            if st.button("Calculate Stoichiometry"):
                known_coeff = reaction.reac.get(known_substance) or reaction.prod.get(known_substance)
                desired_coeff = reaction.reac.get(desired_substance) or reaction.prod.get(desired_substance)

                if known_coeff is None or desired_coeff is None:
                    st.error("Substances not found in the equation.")
                else:
                    # Use mole ratio
                    result = (desired_coeff / known_coeff) * known_moles
                    st.success(f"{result:.4f} mol of {desired_substance} will be produced/consumed.")

                col1, col2 = st.columns(2)
                with col1:
                    st.info(f"**{known_substance}**: {known_moles:.4f} mol")
                    known_substance_mol = df[df["Formula"] == known_substance]["mol"]
                    if not known_substance_mol.empty:
                        st.image(Draw.MolToImage(known_substance_mol.iloc[0]), caption=f"Structure of {known_substance}", use_container_width=True)
                    else:
                        st.warning(f"No structure available for {known_substance}.")
                with col2:
                    st.info(f"**{desired_substance}**: {result:.4f} mol")
                    desired_substance_mol = df[df["Formula"] == desired_substance]["mol"]
                    if not desired_substance_mol.empty:
                        st.image(Draw.MolToImage(desired_substance_mol.iloc[0]), caption=f"Structure of {desired_substance}", use_container_width=True)
                    else:
                        st.warning(f"No structure available for {desired_substance}.")                

    except Exception as e:
        st.error(f"Error parsing equation: {e}")

        # if st.button("Reset"):
        #     st.session_state.show_stoich = False

    

def molarity_calculator():
    with st.expander("📘 What is Molarity?"):
        st.write("""
        Molarity is a measure of concentration defined as the number of moles of solute per liter of solution.
        It is commonly used in chemistry to express concentrations of solutions.
        Molarity = (moles of solute) / (volume of solution in L)
                """)
    moles = st.number_input("Enter moles of solute (mol)", min_value=0.0, format="%.4f")
    volume = st.number_input("Enter volume of solution (L)", min_value=0.0, format="%.4f")
    
    if st.button("Calculate Molarity"):
        if volume > 0:
            molarity = moles / volume
            st.success(f"Molarity = {molarity:.4f} mol/L")
        else:
            st.warning("Volume must be greater than 0.")

    

def molality_calculator():
    with st.expander("📘 What is Molality?"):
        st.write("""
        Molality is a measure of concentration defined as the number of moles of solute per kilogram of solvent.
        It is used in colligative properties and is independent of temperature.
        Molality = (moles of solute) / (mass of solvent in kg)
        """)
    moles = st.number_input("Enter moles of solute (mol)", min_value=0.0, format="%.4f", key="mol_molality")
    mass_solvent = st.number_input("Enter mass of solvent (kg)", min_value=0.0, format="%.4f")

    if st.button("Calculate Molality"):
        if mass_solvent > 0:
            molality = moles / mass_solvent
            st.success(f"Molality = {molality:.4f} mol/kg")
        else:
            st.warning("Solvent mass must be greater than 0.")
    
    

def normality_calculator():

    with st.expander("📘 What is Normality?"):
                st.write("""
                Normality is a measure of concentration equivalent to molarity multiplied by the number of equivalents.
                It is used in acid-base reactions and redox reactions.
                Normality = (gram equivalents of solute) / (volume of solution in L)
                """)

    equivalents = st.number_input("Enter gram equivalents of solute", min_value=0.0, format="%.4f")
    volume = st.number_input("Enter volume of solution (L)", min_value=0.0, format="%.4f", key="vol_normality")

    if st.button("Calculate Normality"):
        if volume > 0:
            normality = equivalents / volume
            st.success(f"Normality = {normality:.4f} eq/L")
        else:
            st.warning("Volume must be greater than 0.")
    
    

def ppm_calculator():
    with st.expander("📘 What is PPM?"):
                st.write("""
                PPM is parts per million, commonly used for very dilute solutions.
                It indicates the mass of solute per million parts of solution.""")

    solute_mass = st.number_input("Mass of solute (mg)", min_value=0.0, format="%.4f")
    solution_volume = st.number_input("Volume of solution (L)", min_value=0.0, format="%.4f", key="vol_ppm")


    if st.button("Calculate PPM"):
        if solution_volume > 0:
            ppm = solute_mass / solution_volume
            st.success(f"PPM = {ppm:.2f} mg/L")
        else:
            st.warning("Volume must be greater than 0.")


    

if __name__ == "__main__":
    display_basic_calculator()