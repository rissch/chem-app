import streamlit as st
import random
import re
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go

from other_pages.edu_tools import generate_equation


st.set_page_config(page_title="Chemistry Educational Tools", layout="wide")
st.title("📘 Chemistry Educational Tools")

# Sidebar Navigation
tool = st.sidebar.selectbox("Select Tool", [
    "Interactive Periodic Table",
    "Chemical Equation Practice",
    "Visual Molecule Builder (Demo)",
    "Concept Visualizer",
    "Flashcards / Quiz"
])

# ---------------- Interactive Periodic Table ----------------
if tool == "Interactive Periodic Table":
    st.subheader("🔬 Interactive Periodic Table")
    st.markdown("Hover over elements to see properties and click for common compounds.")

    # Embed the Ptable interactive periodic table using iframe
    st.subheader("Interactive Periodic Table")
    iframe_html = """
    <iframe src="https://ptable.com/#Properties" width="100%" height="600px"></iframe>
    """
    st.markdown(iframe_html, unsafe_allow_html=True)

    # Load data
    @st.cache_data
    def load_data():
        url = "https://raw.githubusercontent.com/Bowserinator/Periodic-Table-JSON/master/PeriodicTableJSON.json"
        df = pd.read_json(url)
        elements = pd.json_normalize(df['elements'])
        return elements

    st.title("🧪 Interactive Periodic Table")

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




# ---------------- Chemical Equation Practice ----------------


if tool == "Chemical Equation Practice":
    st.subheader("🧪 Practice Balancing Equations")
    equation = generate_equation()
    st.markdown(f"Balance the following reaction:")
    st.code(equation)
    user_input = st.text_input("Your Balanced Equation")
    if user_input:
        st.success("(Demo) Thanks! Check with a chemistry tool or teacher to verify.")

# ---------------- Visual Molecule Builder (Demo) ----------------
if tool == "Visual Molecule Builder (Demo)":
    st.subheader("🧱 Build a Molecule (Demo)")
    atoms = ["C", "H", "O", "N", "Cl"]
    atom = st.selectbox("Choose Atom", atoms)
    count = st.number_input("How many?", min_value=1, max_value=10, step=1)
    structure = atom * count
    st.markdown(f"Molecule SMILES: `{structure}`")
    st.info("(Full visual builder needs integration with RDKit drawing tools)")

# ---------------- Concept Visualizer ----------------
if tool == "Concept Visualizer":
    st.subheader("📊 Chemistry Concept Visualizer")
    concept = st.selectbox("Select Concept", ["Le Chatelier's Principle", "Redox Reactions", "Entropy", "Acid-Base Neutralization"])
    if concept == "Le Chatelier's Principle":
        st.image("le_chatelier.png", caption="Le Chatelier's Principle")
    elif concept == "Redox Reactions":
        st.image("redox_reaction.png", caption="Redox Reaction Example")
    elif concept == "Entropy":
        st.image("entropy_diagram.png", caption="Entropy Increase")
    elif concept == "Acid-Base Neutralization":
        st.image("neutralization.png", caption="Acid-Base Reaction")

# ---------------- Flashcards / Quiz ----------------
flashcards = {
    "What is the molar mass of H2O?": "18.02 g/mol",
    "Name the compound with formula NaCl": "Sodium chloride",
    "What is the symbol for Potassium?": "K",
    "What gas is produced when acid reacts with metal?": "Hydrogen"
}

if tool == "Flashcards / Quiz":
    st.subheader("🧠 Flashcards / Quiz")
    question = random.choice(list(flashcards.keys()))
    st.markdown(f"**Question:** {question}")
    answer = st.text_input("Your Answer")
    if answer:
        correct = flashcards[question].lower()
        if answer.lower() == correct:
            st.success("Correct!")
        else:
            st.error(f"Incorrect. Correct answer: {flashcards[question]}")
