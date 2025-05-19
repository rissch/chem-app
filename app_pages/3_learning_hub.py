import streamlit as st
from streamlit_option_menu import option_menu

def display_learning_hub():
    selected = option_menu(
        menu_title="Learning Hub",
        options=[
            "Chemistry Concepts Tutorials",
            "Problem Solving Guides",
            "Downloadable Reference Materials",
            "Visuals and Tutorials",
        ],
        icons=[
            "book", "calculator", "file-earmark-text", "play-circle"
        ],
        menu_icon="book",
        default_index=0,
        orientation="horizontal",
    )

    if selected == "Chemistry Concepts Tutorials":
        # Page title
        st.title("📘 Chemistry Concept Tutorials")

        st.markdown("Welcome! Here's a quick and simple guide to help you understand some of the most important concepts in chemistry.")

        # Dropdown menu for selecting topic
        concept = st.selectbox(
            "📚 Select a Topic to Learn",
            [
                "Stoichiometry",
                "Chemical Equilibrium",
                "Reaction Kinetics",
                "Acids and Bases"
            ]
        )

        # Stoichiometry
        if concept == "Stoichiometry":
            st.header("⚖️ Stoichiometry Made Easy")
            st.markdown("""
            **Stoichiometry** is the calculation of **reactants and products** in a chemical reaction.

            **Why it's important:** It helps you predict how much product you'll get or how much reactant you need.

            ### Example:
            **Reaction:** `2H₂ + O₂ → 2H₂O`

            This means:
            - 2 moles of hydrogen react with 1 mole of oxygen to make 2 moles of water.

            ### Basic Steps:
            1. **Balance the equation.**
            2. **Convert grams → moles** (use molar mass).
            3. **Use mole ratio** from the balanced equation.
            4. **Convert moles → grams** (if needed).

            ✅ **Tip:** Always double-check units and mole ratios!
            """)

            # Chemical Equilibrium
        elif concept == "Chemical Equilibrium":
                st.header("🔁 Chemical Equilibrium Explained")
                st.markdown("""
            **Equilibrium** is when a reversible reaction’s forward and backward rates become **equal**.

            **Example:**  
            `N₂ + 3H₂ ⇌ 2NH₃`  
            After some time, the amount of ammonia (NH₃) stops increasing — not because the reaction stopped, but because it's **going both ways at the same speed**.

            ### Key Concepts:
            - **Dynamic**: Molecules still react.
            - **Constant concentrations**: But no net change.
            - **Le Chatelier’s Principle**: If conditions change (like pressure or temperature), the reaction will shift to restore balance.

            ✅ **Tip:** Think of equilibrium like a tug-of-war with equal strength on both sides.
            """)

            # Reaction Kinetics
        elif concept == "Reaction Kinetics":
            st.header("⏱️ Understanding Reaction Kinetics")
            st.markdown("""
            **Kinetics** is the study of how **fast** reactions happen.

            ### Factors That Affect Speed:
            1. **Concentration** – More particles = more collisions.
            2. **Temperature** – Hotter = faster particles = more energy.
            3. **Surface Area** – More exposed area = more collisions.
            4. **Catalysts** – Speed up reaction without being used up.

            ### Rate of Reaction Formula:
            `Rate = change in concentration / time`

            ✅ **Tip:** Faster reactions usually mean more frequent and energetic collisions!
            """)

            # Acids and Bases
        elif concept == "Acids and Bases":
            st.header("🧪 Acids & Bases Demystified")
            st.markdown("""
            ### Definitions:
            - **Acid**: Donates H⁺ (proton)
            - **Base**: Accepts H⁺

            ### Common Examples:
            - Acid: HCl, H₂SO₄
            - Base: NaOH, NH₃

            ### pH Scale:
            - **0–6** = Acidic  
            - **7** = Neutral  
            - **8–14** = Basic

            ### Neutralization Reaction:
            `Acid + Base → Salt + Water`

            ### Strong vs. Weak:
            - **Strong acids/bases** fully dissociate in water.
            - **Weak** ones partially dissociate.

            ✅ **Tip:** Use indicators (like litmus or phenolphthalein) to test for acidity/basicity!
            """)

            # Closing message
            st.markdown("---")
            st.success("Want more topics like thermochemistry, redox, or atomic structure? Let me know!")
    elif selected == "Problem Solving Guides":
        st.title("🧠 Chemistry Problem Solving Guides")
        st.markdown("Master problem-solving in chemistry with easy-to-follow strategies for every topic!")

        # Guide selection
        problem_type = st.selectbox("📘 Choose a Problem Type:", [
            "Mole Concept",
            "Limiting Reagent",
            "Empirical & Molecular Formula",
            "Concentration (Molarity)",
            "pH Calculations",
            "Gas Laws (Ideal Gas Equation)",
            "Thermochemistry",
        ])

        # Guide Content
        if problem_type == "Mole Concept":
            st.header("⚖️ Mole Concept Problem Guide")
            st.markdown("""
        **Goal:** Convert between grams, moles, and number of particles.

        ### 🛠️ Steps:
        1. **Know the molar mass** (g/mol).
        2. Use formulas:
        - `Moles = mass / molar mass`
        - `Particles = moles × Avogadro's number (6.022×10²³)`

        ### ✅ Example:
        **Q:** How many moles are in 36g of water (H₂O)?  
        **A:** Molar mass = 18 g/mol → Moles = 36 / 18 = **2 mol**
        """)

        elif problem_type == "Limiting Reagent":
            st.header("🧪 Limiting Reagent Problem Guide")
            st.markdown("""
        **Goal:** Find which reactant gets used up first in a reaction.

        ### 🛠️ Steps:
        1. **Balance the chemical equation.**
        2. Convert all reactant amounts to moles.
        3. Divide moles by their coefficients.
        4. The **smallest result** is the limiting reagent.

        ### ✅ Example:
        Reacting 4 mol H₂ and 2 mol O₂  
        Equation: `2H₂ + O₂ → 2H₂O`  
        → H₂: 4/2 = 2, O₂: 2/1 = 2 → **Limiting: Neither** (perfect ratio)
        """)

        elif problem_type == "Empirical & Molecular Formula":
            st.header("📊 Empirical & Molecular Formula Guide")
            st.markdown("""
        **Empirical Formula:** Simplest whole-number ratio of atoms.  
        **Molecular Formula:** Actual number of atoms in a molecule.

        ### 🛠️ Steps:
        1. Convert mass % to grams (assume 100g).
        2. Convert grams to moles.
        3. Divide all moles by the smallest.
        4. Multiply to get whole numbers (if needed).

        ### ✅ Example:
        70% Fe, 30% O  
        → 70g Fe / 55.85 ≈ 1.25 mol  
        → 30g O / 16 ≈ 1.875 mol  
        → Ratio Fe:O = 1.25 : 1.875 → Divide both by 1.25 → 1 : 1.5 → Multiply by 2 → **Fe₂O₃**
        """)

        elif problem_type == "Concentration (Molarity)":
            st.header("🧫 Molarity Problem Guide")
            st.markdown("""
        **Molarity (M)** = moles of solute / liters of solution

        ### 🛠️ Steps:
        1. Convert mass to moles if needed.
        2. Convert mL to liters.
        3. Use `M = mol / L`

        ### ✅ Example:
        Q: What’s the molarity of 0.5 mol NaCl in 250 mL solution?  
        A: 250 mL = 0.250 L → M = 0.5 / 0.250 = **2.0 M**
        """)

        elif problem_type == "pH Calculations":
            st.header("🧪 pH & pOH Problem Guide")
            st.markdown("""
        ### 🧮 Formulas:
        - `pH = -log[H⁺]`
        - `pOH = -log[OH⁻]`
        - `pH + pOH = 14`

        ### 🛠️ Steps:
        1. Use concentration of H⁺ or OH⁻.
        2. Apply the log formula.
        3. Use the relationship if needed.

        ### ✅ Example:
        Q: What is the pH of 0.01 M HCl?  
        A: pH = -log(0.01) = 2
        """)

        elif problem_type == "Gas Laws (Ideal Gas Equation)":
            st.header("🌬️ Ideal Gas Equation Guide")
            st.markdown("""
        **Equation:** PV = nRT

        Where:
        - P = pressure (atm)
        - V = volume (L)
        - n = moles
        - R = 0.0821 L·atm/mol·K
        - T = temperature (K)

        ### 🛠️ Steps:
        1. Convert temperature to Kelvin (K = °C + 273).
        2. Rearrange formula to solve for unknown.

        ### ✅ Example:
        Q: What is the volume of 2 mol gas at 1 atm and 273 K?  
        A: V = nRT/P = (2 × 0.0821 × 273) / 1 ≈ **44.8 L**
        """)

        elif problem_type == "Thermochemistry":
            st.header("🔥 Thermochemistry Guide")
            st.markdown("""
        **Goal:** Understand heat changes in reactions.

        ### 🛠️ Key Formulas:
        - `q = mcΔT`
        - `ΔH = q / mol`

        ### Definitions:
        - **q** = heat (J)
        - **m** = mass (g)
        - **c** = specific heat (4.18 J/g·°C for water)
        - **ΔT** = change in temperature

        ### ✅ Example:
        Q: Heat absorbed by 100g water heated from 25°C to 75°C?  
        A: q = 100 × 4.18 × (75 - 25) = **20,900 J**
        """)

        # Footer
        st.markdown("---")
        st.success("Keep practicing! Chemistry becomes easier with clear steps and consistent problem solving.")
    elif selected == "Downloadable Reference Materials":
        st.title("📥 Downloadable Reference Materials")

        if st.button("Download Basic Chemistry Tutorials & Visuals"):
            with open("downloadable_reference_files/Basic Chemistry Tutorials & Visuals.pdf", "rb") as file:
                st.info("Preparing file for download...")
                file_data = file.read()
            
            st.success("File ready for download!")
            st.download_button(
                label="Download Example PDF",
                data=file_data,
                file_name="Basic Chemistry Tutorials & Visuals.pdf",
                mime="application/pdf",
            )

        st.markdown("---")

        with open("downloadable_reference_files/Exam Preparation tips.pdf", "rb") as file:
            file_data = file.read()
        st.download_button(
            label="Download Exam Preparation tips",
            data=file_data,
            file_name="Exam Preparation tips.pdf",
            mime="application/pdf",
        )
    elif selected == "Visuals and Tutorials":
        st.title("🔬 Basic Chemistry Tutorials & Visuals")
        st.markdown("""
        Welcome to the Basic Chemistry learning hub!  
        Explore videos, animations, and visuals for core chemistry topics to boost your understanding.
        """)

        # --- Section: Stoichiometry ---
        st.header("1. Stoichiometry")
        st.markdown("""
        **What is Stoichiometry?**  
        Stoichiometry deals with the quantitative relationships between reactants and products in chemical reactions.
        """)

        st.video("https://www.youtube.com/watch?v=Gle1bPAZsgg")  # Example Stoichiometry video

        # st.image(
        #     "https://upload.wikimedia.org/wikipedia/commons/thumb/7/72/Reaction_Stoichiometry.svg/800px-Reaction_Stoichiometry.svg.png",
        #     caption="Stoichiometry Concept Diagram",
        #     use_container_width=True
        # )

        st.markdown("### 🔬 Stoichiometry Animation")
        st.markdown(
            "[Click here to view the animation](https://www.dlt.ncssm.edu/core/Chapter6-Stoichiometry/Chapter6-Animations/Stoichiometry-1.html)",
            unsafe_allow_html=True
        )

        st.markdown("---")

        # --- Section: Chemical Equilibrium ---
        st.header("2. Chemical Equilibrium")
        st.markdown("""
        **Chemical Equilibrium** occurs when the rates of forward and reverse reactions are equal, resulting in stable concentrations of reactants and products.
        """)

        st.video("https://www.youtube.com/watch?v=JsoawKguU6A")  # Example Equilibrium video

        # # Embed interactive simulation from PhET (iframe)
        # equilibrium_iframe = """
        # <div style="width:700px; height:400px;">
        # <iframe src="https://phet.colorado.edu/sims/html/equilibrium/latest/equilibrium_en.html" width="700" height="400" style="border:none;"></iframe>
        # </div>
        # """
        

        # components.html(equilibrium_iframe, height=420)

        st.markdown("---")

        # --- Section: Reaction Kinetics ---
        st.header("3. Reaction Kinetics")
        st.markdown("""
        **Reaction Kinetics** studies the speed or rate of chemical reactions and the factors affecting them.
        """)

        st.video("https://www.youtube.com/watch?v=DbEecMtO6aY")  # Example Kinetics video

        # st.image(
        #     "https://chem.libretexts.org/@api/deki/files/4372/Kinetics_graph.png?revision=1",
        #     caption="Reaction Rate Graph",
        #     use_container_width=True
        # )

        st.markdown("""
        <div style="background-color: white; padding: 30px; border-radius: 10px; text-align: center;">
            <img src="https://chem.libretexts.org/@api/deki/files/4372/Kinetics_graph.png?revision=1", 
                >
        </div>
        """, unsafe_allow_html=True)

        st.markdown("---")

        # --- Section: Acids and Bases ---
        st.header("4. Acids and Bases")
        st.markdown("""
        Acids and bases are substances that donate or accept protons (H⁺ ions), respectively, and are fundamental to many chemical processes.
        """)

        st.video("https://www.youtube.com/watch?v=V5Mq_cL9Bck")  # Example Acids and Bases video

        st.image(
            "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0b/Ph_scale.svg/800px-Ph_scale.svg.png",
            caption="pH Scale Chart",
            use_container_width=True
        )

        st.markdown("---")


if __name__ == "__main__":
    display_learning_hub()