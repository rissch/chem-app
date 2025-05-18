from collections import defaultdict
import math
import random
import re
import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Draw
from streamlit.components.v1 import html
import json

from utils import get_atoms


def display_edu_tools():

    html("""
<script>
    const screenWidth = window.innerWidth;
    console.log("Screen width detected:", screenWidth);
    const params = new URLSearchParams(window.location.search);
    if (!params.get("width")) {
        params.set("width", screenWidth);
        console.log("Reloading with width param:", params.toString());
        window.location.search = params.toString();
    } else {
        console.log("Width param already present:", params.get("width"));
    }
    </script>
    """, height=0)

    # ✅ Step 2: Access query parameters correctly
    query_params = st.query_params  # note: no parentheses

    # ✅ Fallback if 'width' not present
    width = int(query_params.get("width", "1000"))

    # ✅ Set orientation based on width
    orientation = "vertical" if width < 768 else "horizontal"
    calc = option_menu(
        menu_title="Educational Tools",
        options=["Periodic Table", "Chemical Equations Practice", "Visual Molecule Builder", "Flashcards"],
        icons=["table", "book", "lightbulb", "question-circle"],
        orientation=orientation
    )
    if calc == "Periodic Table":
        st.title("🧪 Interactive Periodic Table")
        st.markdown("Hover over elements to see properties and click for common compounds.")

        # Embed the Ptable interactive periodic table using iframe
        # iframe_html = """
        # <iframe src="https://ptable.com/#Properties" width="100%" height="600px"></iframe>
        # """
        # st.markdown(iframe_html, unsafe_allow_html=True)

        # Load data
        @st.cache_data
        def load_data():
            file_path = "PeriodicTableJSON.json"
            df = pd.read_json(file_path)
            elements = pd.json_normalize(df['elements'])
            return elements


        df = load_data()

        # Positioning
        df['x'] = df['xpos']
        df['y'] = -df['ypos']

        # Better color map
        category_colors = {
            'diatomic nonmetal': '#1f77b4',
            'noble gas': '#9467bd',
            'alkali metal': '#ff7f0e',
            'alkaline earth metal': '#2ca02c',
            'metalloid': '#d62728',
            'polyatomic nonmetal': '#17becf',
            'transition metal': '#8c564b',
            'post-transition metal': '#bcbd22',
            'lanthanide': '#e377c2',
            'actinide': '#7f7f7f',
            'unknown, probably transition metal': '#008080', 
            'unknown, probably metalloid': '#FF7F50',  
            'unknown, probably post-transition metal': '#5A3332',  
            'unknown, predicted to be noble gas': '#E2805D',  
            'unknown, but predicted to be an alkali metal': '#77702A',  
            '': '#cccccc'
        }

        df['color'] = df['category'].map(category_colors).fillna("#dddddd")

        # Plot
        fig = go.Figure()

        for _, row in df.iterrows():
            # Adjust text color for contrast
            text_color = 'white' if row['category'] not in ['unknown', 'post-transition metal'] else 'black'
            fig.add_trace(go.Scatter(
                x=[row['x']],
                y=[row['y']],
                mode='markers+text',
                marker=dict(size=42, color=row['color'], line=dict(width=1, color='black')),  # Smaller circles
                text=row['symbol'],
                textfont=dict(size=14, color=text_color),
                hovertemplate=f"<b>{row['name']}</b><br>Atomic #: {row['number']}<br>Mass: {row['atomic_mass']}<br>Category: {row['category']}<extra></extra>"
            ))

        # Layout
        fig.update_layout(
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            showlegend=False,
            height=800,
            margin=dict(t=20, b=20, l=20, r=20),
            plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
            paper_bgcolor='rgba(0,0,0,0)'  # Transparent canvas
        )


        st.plotly_chart(fig, use_container_width=True)

    

    # In your Streamlit display function:
    elif calc == "Chemical Equations Practice":
        # function defination
        with open("ChemicalEquations.json") as f:
            eq_dict = json.load(f)
        
        def generate_equation():
        # A list of unbalanced equations (demo purposes)
            return eq_dict[0]['balance_reaction_examples']

        def parse_compound(compound):
            elements = re.findall(r'([A-Z][a-z]?)(\d*)', compound)
            count_dict = defaultdict(int)
            for element, count in elements:
                count_dict[element] += int(count) if count else 1
            return count_dict

        def parse_side(side):
            terms = side.split('+')
            total_counts = defaultdict(int)
            for term in terms:
                term = term.strip()
                match = re.match(r"(\d*)\s*([A-Za-z0-9]+)", term)
                if not match:
                    continue
                coeff = int(match.group(1)) if match.group(1) else 1
                compound = match.group(2)
                atom_counts = parse_compound(compound)
                for element, count in atom_counts.items():
                    total_counts[element] += coeff * count
            return total_counts

        def is_balanced(user_eq):
            try:
                if '->' not in user_eq:
                    return False
                left, right = user_eq.split('->')
                left_counts = parse_side(left)
                right_counts = parse_side(right)
                return left_counts == right_counts
            except Exception:
                return False
        
        
        
        st.subheader("🧪 Practice Balancing Equations")

        with st.sidebar:
            selected_option = st.radio("Select Option : ", ["Random", "Custom"], key="eq_option")

        eq_dict = generate_equation()

        # -- Session state for random equation --
        if "random_eq" not in st.session_state:
            st.session_state.random_eq = random.choice(list(eq_dict.items()))

        # -- Button to shuffle the equation only in random mode --
        if selected_option == "Random":
            if st.button("🔀 Shuffle Equation"):
                st.session_state.random_eq = random.choice(list(eq_dict.items()))

        # -- Equation selection logic --
        if selected_option == "Custom":
            selected_equation = st.selectbox("Select an equation to balance", list(eq_dict.keys()), key="eq_select")
            unbalanced = selected_equation
            balanced = eq_dict[selected_equation]
        else:
            unbalanced, balanced = st.session_state.random_eq

        # -- Display the equation --
        st.markdown("**Balance the following chemical equation:**")
        st.code(unbalanced)

        # -- User input and validation --
        user_input = st.text_input("Enter your balanced equation")

        if user_input:
            if is_balanced(user_input) and user_input.replace(" ", "") == balanced.replace(" ", ""):
                st.success("✅ Well done! Your equation is balanced correctly.")
            elif is_balanced(user_input):
                st.warning("⚠️ That equation is balanced, but it's not the correct one selected reaction.")
            else:
                st.error("❌ That doesn't seem correct. Try again or check your syntax.")
            
            with st.expander("See the correct balanced equation"):
                st.code(balanced)

    elif calc == "Visual Molecule Builder":

        st.subheader("🧱 Visual Molecule Builder")

        atoms = get_atoms()
        atom = st.selectbox("Choose Atom", atoms)
        atom = atom[0]  # Get the symbol of the selected atom
        count = st.number_input("How many atoms?", min_value=1, max_value=6, step=1)


        # Build a linear molecule with single bonds (e.g., C-C-C-C...)
        mol = Chem.RWMol()

        atom_indices = []
        for _ in range(count):
            idx = mol.AddAtom(Chem.Atom(atom))
            atom_indices.append(idx)

        # Connect atoms with single bonds
        for i in range(count - 1):
            mol.AddBond(atom_indices[i], atom_indices[i + 1], Chem.BondType.SINGLE)

        final_mol = mol.GetMol()
        smiles = Chem.MolToSmiles(final_mol)
        img = Draw.MolToImage(final_mol, size=(300, 300))

        st.image(img, caption=f"Molecule Structure")
        st.markdown(f"**SMILES String**: `{smiles}`")


    # elif calc == "Concept Visualizer":
    #     st.title("🧠 Concept Visualizer")

    #     # 1. Atomic Structure
    #     with st.expander("1. 🧬 Atomic Structure"):
    #         st.markdown("Visualize a simple Bohr model")
    #         atomic_number = st.slider("Atomic Number", 1, 20, 1)
    #         shells = [2, 8, 8]  # simplified configuration
    #         electrons = atomic_number
    #         layout = []
    #         for i, max_e in enumerate(shells):
    #             radius = (i + 1) * 20
    #             count = min(electrons, max_e)
    #             for j in range(count):
    #                 angle = 2 * math.pi * j / count
    #                 x = radius * math.cos(angle)
    #                 y = radius * math.sin(angle)
    #                 layout.append(go.Scatter(x=[x], y=[y], mode='markers', marker=dict(size=10), showlegend=False))
    #             electrons -= count
    #         fig = go.Figure(data=layout)
    #         fig.update_layout(height=400, xaxis=dict(visible=False), yaxis=dict(visible=False), plot_bgcolor='white')
    #         st.plotly_chart(fig, use_container_width=True)

        # # 2. Electron Configuration
        # with st.expander("2. ⚛️ Electron Configuration"):
        #     st.markdown("Show orbital filling order for atomic number")
        #     configs = ["1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p"]
        #     z = st.slider("Atomic Number", 1, 30, 1)
        #     fill_order = []
        #     max_electrons = {"s": 2, "p": 6, "d": 10}
        #     electrons = z
        #     for orb in configs:
        #         subshell = orb[-1]
        #         max_e = max_electrons[subshell]
        #         count = min(electrons, max_e)
        #         if count > 0:
        #             fill_order.append(f"{orb}^{count}")
        #         electrons -= count
        #         if electrons <= 0:
        #             break
        #     st.write("Electron Configuration:", " ".join(fill_order))

        # # 3. Periodic Trends
        # with st.expander("3. 📊 Periodic Trends"):
        #     st.markdown("Visualize atomic radius across period 2")
        #     elements = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
        #     radii = [167, 112, 87, 67, 56, 48, 42, 38]  # pm
        #     fig = go.Figure([go.Bar(x=elements, y=radii)])
        #     fig.update_layout(title="Atomic Radius (pm)", yaxis_title="Radius (pm)")
        #     st.plotly_chart(fig, use_container_width=True)

        # # 4. Molecular Orbitals
        # with st.expander("4. 🌀 Molecular Orbitals"):
        #     st.markdown("Demo only — actual visuals require SVGs")
        #     st.image("https://chem.libretexts.org/@api/deki/files/18388/clipboard_e3a49c68c9fc3df5ef8c27277fd8c3fd1.png")

        # # 5. Bonding Types
        # with st.expander("5. 🧲 Bonding Types"):
        #     st.markdown("Textual comparison with examples")
        #     st.write("**Ionic:** NaCl\n\n**Covalent:** H2O\n\n**Metallic:** Cu")

        # # 6. Hybridization
        # with st.expander("6. 🧪 Hybridization"):
        #     st.markdown("Visual of sp3 tetrahedral molecule")
        #     st.image("https://chem.libretexts.org/@api/deki/files/21824/sp3_Hybridization.png")

        # # 7. States of Matter
        # with st.expander("7. 🌡️ States of Matter"):
        #     st.markdown("Dynamic view of particles in different states")
        #     st.image("https://cdn.britannica.com/88/20688-050-7A01C59B/states-of-matter.jpg")

        # # 8. Reaction Energy Diagrams
        # with st.expander("8. 🔥 Reaction Energy Diagrams"):
        #     st.markdown("Exothermic reaction profile")
        #     fig = go.Figure()
        #     fig.add_trace(go.Scatter(x=[0, 1, 2], y=[0, 5, -3], mode='lines+markers'))
        #     fig.update_layout(title="Exothermic Reaction", xaxis_title="Reaction Progress", yaxis_title="Energy")
        #     st.plotly_chart(fig, use_container_width=True)

        # # 9. Le Chatelier’s Principle
        # with st.expander("9. ⚖️ Le Chatelier’s Principle"):
        #     st.markdown("Interactive shift with concentration change")
        #     conc = st.slider("Increase Reactant Concentration", 1.0, 5.0, 1.0)
        #     st.write(f"Equilibrium shifts to products by factor {conc}")

        # # 10. Acid-Base Concepts
        # with st.expander("10. 🧪 Acid-Base Concepts"):
        #     st.markdown("pH Scale Visual")
        #     fig = go.Figure([go.Bar(x=list(range(1, 15)), y=[1]*14, marker_color=list(range(1, 15)))])
        #     fig.update_layout(title="pH Scale", xaxis_title="pH", yaxis=dict(showticklabels=False))
        #     st.plotly_chart(fig, use_container_width=True)

    elif calc == "Flashcards":

        # Flashcard data structure
        flashcards = {
            "Atomic Structure": {
                "question": "What are the three subatomic particles?",
                "answer": "Proton (+), Neutron (0), Electron (-)",
                "details": "Protons and neutrons reside in the nucleus, while electrons orbit around it. Atomic number = protons; mass number = protons + neutrons."
            },
            "Periodic Table": {
                "question": "What determines an element’s position in the periodic table?",
                "answer": "Atomic number",
                "details": "Elements are arranged in order of increasing atomic number, which represents the number of protons in the nucleus."
            },
            "Chemical Bonding": {
                "question": "What is the difference between ionic and covalent bonds?",
                "answer": "Ionic: Transfer of electrons; Covalent: Sharing of electrons",
                "details": "Ionic bonds occur between metals and non-metals, covalent between non-metals."
            },
            "Moles and Molar Mass": {
                "question": "How many particles are in one mole?",
                "answer": "6.022 × 10²³ particles",
                "details": "Avogadro’s number represents the number of atoms, molecules, or ions in a mole."
            },
            "Stoichiometry": {
                "question": "What does stoichiometry help calculate?",
                "answer": "Quantities of reactants and products",
                "details": "It uses balanced chemical equations to determine the ratios and amounts."
            },
            "Acids and Bases": {
                "question": "What is the pH of a neutral solution?",
                "answer": "7",
                "details": "pH < 7: Acidic; pH = 7: Neutral; pH > 7: Basic"
            },
            "Gas Laws": {
                "question": "What law states that volume is inversely proportional to pressure?",
                "answer": "Boyle’s Law",
                "details": "P1V1 = P2V2 at constant temperature"
            },
            "Thermodynamics": {
                "question": "What is the first law of thermodynamics?",
                "answer": "Energy cannot be created or destroyed.",
                "details": "ΔU = q + w, where U is internal energy, q is heat, and w is work."
            },
            "Kinetics": {
                "question": "What does reaction rate depend on?",
                "answer": "Concentration, temperature, surface area, and catalysts",
                "details": "Higher concentration or temperature usually increases the reaction rate."
            },
            "Equilibrium": {
                "question": "What happens to a system at equilibrium?",
                "answer": "Rate of forward and reverse reactions is equal",
                "details": "Concentrations of reactants and products remain constant over time."
            }
        }

        st.title("📚 Chemistry Flashcards")

        # Flashcard selector
        concept = st.selectbox("Choose a Concept", list(flashcards.keys()))

        # Display the flashcard
        card = flashcards[concept]
        with st.expander(f"❓ {card['question']}"):
            st.write(f"✅ **Answer:** {card['answer']}")
            if st.toggle("Show Detailed Explanation"):
                st.markdown(f"📝 {card['details']}")


elist = ["H2", "O2", "N2", "Cl2", "Na", "K", "Mg", "Ca", "Al", "Fe", "Cu", "Zn", "HCl", "NaOH", "H2SO4"]
def generate_equation():
    reac = random.sample(elist, 2)
    prod = random.sample(elist, 1) + ["H2O"]
    return f"{reac[0]} + {reac[1]} -> {prod[0]} + {prod[1]}"