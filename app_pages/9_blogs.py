import streamlit as st
from streamlit_option_menu import option_menu

def display_blogs():
    with st.sidebar:
        selected_blog_option = st.radio(
            "Blog Categories",
            options=[
                "General Blogs",
                "On Research Blogs",],
            index=0
        )
    
    if selected_blog_option == "General Blogs":
        selected = option_menu(
            menu_title="Blogs",
            options=[
                "Study Tips",
                "Exam Prepration",
                "Lab Safety",
                "Career Guides"
                ],
                orientation="horizontal",
                icons=[
                    "book", "check-circle", "shield-lock", "briefcase"
                ],
                default_index=0,)
                
        
        if selected == "Study Tips":
            st.title("🧪 Mastering Chemistry: Study Tips for Success")

            st.markdown("""
            Chemistry is often called the *central science* — it connects biology, physics, environmental science, and even engineering. But for many students, it can feel like a maze of confusing formulas, abstract concepts, and endless reactions.

            Whether you're a high school student, a pre-med hopeful, or just someone fascinated by the building blocks of matter, here are **practical and powerful study tips** to help you conquer chemistry with confidence.
            """)

            # Tip 1
            st.subheader("🔬 1. Understand the Why, Not Just the What")
            st.markdown("""
            Chemistry isn't just about memorizing reactions — it's about understanding *why* atoms behave the way they do.

            **Ask questions like:**
            - Why does sodium react violently with water?
            - Why is carbon the backbone of organic chemistry?
            - Why does increasing temperature speed up a reaction?

            ✅ **Tip:** Whenever you learn a concept, link it to real-life examples — this will boost both retention and curiosity.
            """)

            # Tip 2
            st.subheader("🧠 2. Master the Basics First")
            st.markdown("""
            A strong foundation makes everything easier. Make sure you're confident with:
            - Periodic trends (groups, periods)
            - Atomic structure (orbitals, electrons)
            - Bonding (ionic, covalent, metallic)
            - Moles, molarity, and stoichiometry

            ✅ **Tip:** Use tools like our **Basic Calculator** and **Unit Converter** to reinforce fundamentals.
            """)

            # Tip 3
            st.subheader("📝 3. Practice Problem-Solving Regularly")
            st.markdown("""
            Chemistry is a *doing* science. You can't just read — you need to solve.

            - Work through textbook questions daily
            - Revisit past papers
            - Use active recall and spaced repetition

            ✅ **Tip:** Don’t avoid hard problems — focus on your weakest areas to grow.
            """)

            # Tip 4
            st.subheader("🎨 4. Use Visual Aids and Color Coding")
            st.markdown("""
            Visual learning boosts memory and understanding.

            - Draw mechanisms and reactions
            - Highlight key ideas in color
            - Watch visual explainers like CrashCourse or NileRed

            ✅ **Tip:** Use mind maps to organize concepts like bonding, acids, or kinetics.
            """)

            # Tip 5
            st.subheader("🧪 5. Do Virtual or Real Experiments")
            st.markdown("""
            Seeing is understanding.

            - Try simulations like **PhET**
            - Observe everyday chemistry (rusting, cooking, cleaning)
            - Perform safe home experiments (e.g., vinegar + baking soda)

            ✅ **Tip:** Link every theory to a real-world application.
            """)

            # Tip 6
            st.subheader("🧩 6. Connect Topics Together")
            st.markdown("""
            Don’t treat chemistry like isolated chapters.

            - Thermodynamics explains equilibrium
            - Organic reactions link to acid-base theory
            - Electrochemistry builds on redox

            ✅ **Tip:** At the end of each topic, map how it connects to others.
            """)

            # Tip 7
            st.subheader("🤝 7. Study in Groups (But Wisely)")
            st.markdown("""
            Group study works well when done right:

            - Teach each other tricky topics
            - Solve problems together
            - Quiz one another

            ✅ **Tip:** Teaching is one of the best ways to learn.
            """)

            # Tip 8
            st.subheader("📱 8. Use Technology to Your Advantage")
            st.markdown("""
            Apps and digital tools can supercharge your learning.

            - Khan Academy for video lessons
            - ChemCollective for virtual labs
            - Our **Chemistry Hub App** for calculators and converters

            ✅ **Tip:** Use tech not just for notes, but for interaction and simulation.
            """)

            # Tip 9
            st.subheader("📆 9. Create a Study Schedule")
            st.markdown("""
            Don’t cram — plan.

            - Break topics into chunks
            - Mix problem-solving with revision
            - Use tools like Pomodoro (25m study / 5m break)

            ✅ **Tip:** Consistency beats intensity. Small efforts daily = big results.
            """)

            # Tip 10
            st.subheader("❤️ 10. Fall in Love with Chemistry")
            st.markdown("""
            Make it personal and exciting.

            - Read chemistry stories: why do onions make you cry?
            - Watch documentaries like *Cosmos*
            - Join forums or chemistry clubs

            ✅ **Tip:** Curiosity is the best teacher. Let it drive you forward.
            """)

            # Final Thoughts
            st.markdown("""---""")
            st.header("🎯 Final Thoughts")
            st.markdown("""
            Chemistry may be challenging, but it's one of the most fascinating subjects you'll ever study. It helps you **think logically**, **solve real-world problems**, and **understand the world** at the molecular level.

            With the right strategies, tools, and mindset, you can master it — and even enjoy it.

            Keep experimenting. Keep asking why. And remember: **every small reaction leads to a big change — just like in learning.**
            """)

        elif selected == "Exam Prepration":
            # Page title
            st.title("📚 Smart Exam Preparation: Strategies for Success")

            st.markdown("""
            Exams can be stressful, but with the right strategies, they become an opportunity to **showcase your knowledge** and growth.

            Whether you're preparing for school finals, university exams, or competitive tests, these tips will help you study smarter and perform better.
            """)

            # Tip 1
            st.subheader("🗂 1. Organize Your Syllabus")
            st.markdown("""
            Start by **breaking your syllabus** into smaller, manageable sections. Understand what topics carry the most weight and prioritize accordingly.

            ✅ **Tip:** Create a study map or checklist so you can track your progress and stay motivated.
            """)

            # Tip 2
            st.subheader("📆 2. Create a Realistic Study Schedule")
            st.markdown("""
            Make a daily or weekly schedule that includes:
            - Specific topics to study
            - Time for revision
            - Short breaks

            Avoid cramming. Aim for **consistency over intensity**.

            ✅ **Use Pomodoro technique**: 25 minutes of focused study, 5-minute break.
            """)

            # Tip 3
            st.subheader("🧠 3. Use Active Learning Methods")
            st.markdown("""
            Don't just read — **engage** with the material.

            - Summarize notes in your own words
            - Teach concepts to someone else
            - Use flashcards and mind maps

            ✅ **Tip:** Practice explaining a topic without looking at your notes.
            """)

            # Tip 4
            st.subheader("📝 4. Practice Past Papers")
            st.markdown("""
            Nothing beats exam-style practice.

            - Solve previous years’ papers
            - Time yourself
            - Review your mistakes and weak areas

            ✅ **Tip:** Simulate real exam conditions to reduce anxiety on test day.
            """)

            # Tip 5
            st.subheader("💡 5. Focus on Understanding, Not Memorizing")
            st.markdown("""
            Memorization fades, but understanding lasts.

            - Ask “Why?” and “How?” often
            - Relate topics to real-life scenarios
            - Use diagrams and analogies

            ✅ **Tip:** Build strong foundational understanding before diving into complex topics.
            """)

            # Tip 6
            st.subheader("📴 6. Limit Distractions")
            st.markdown("""
            Turn off notifications, silence your phone, and create a clean study space.

            ✅ **Tip:** Use apps like **Forest** or **Focus To-Do** to stay distraction-free.
            """)

            # Tip 7
            st.subheader("🤝 7. Group Study (if it works for you)")
            st.markdown("""
            Studying with peers can:
            - Clarify difficult concepts
            - Help you discover gaps in your knowledge
            - Keep you motivated

            ✅ **Tip:** Avoid off-topic chatter. Set a specific goal for each session.
            """)

            # Tip 8
            st.subheader("😴 8. Sleep and Eat Well")
            st.markdown("""
            Your brain needs fuel and rest.

            - Sleep at least 7–8 hours
            - Stay hydrated
            - Eat brain foods like nuts, berries, and fish

            ✅ **Tip:** Avoid late-night cramming the day before the exam — it hurts more than it helps.
            """)

            # Tip 9
            st.subheader("🎯 9. Stay Positive and Manage Stress")
            st.markdown("""
            Your mindset matters.

            - Practice deep breathing or light meditation
            - Take breaks when overwhelmed
            - Celebrate small wins

            ✅ **Tip:** Visualize success. Believe in your preparation and effort.
            """)

            # Final Thoughts
            st.markdown("---")
            st.header("🌟 Final Thoughts")
            st.markdown("""
            Exams are **not a measure of your worth**, but they are an opportunity to learn discipline, focus, and perseverance.

            With the right planning, habits, and attitude, **you can perform at your best** and turn stress into success.

            💬 *“Success is the sum of small efforts repeated day in and day out.”* — Robert Collier
            """)

        elif selected == "Lab Safety":
            # Page title
            st.title("🧪 Lab Safety Essentials: Protect Yourself and Others")

            st.markdown("""
            Working in a lab can be exciting and educational, but it also comes with responsibilities. Whether you're a student, researcher, or professional, **lab safety** is essential to prevent accidents and ensure a productive environment.

            This guide will walk you through the **most important lab safety practices** every chemistry student or lab worker must follow.
            """)

            # Section 1
            st.subheader("🔍 1. Know Your Surroundings")
            st.markdown("""
            Before starting any experiment:

            - Familiarize yourself with **emergency exits**, **eye wash stations**, **fire extinguishers**, and **first aid kits**.
            - Understand how to use all lab equipment properly.
            - Keep your workspace clean and organized.

            ✅ **Tip:** Ask for a lab tour if you're new to the space.
            """)

            # Section 2
            st.subheader("👕 2. Wear Proper Lab Attire")
            st.markdown("""
            Always wear:

            - **Lab coat**
            - **Safety goggles**
            - **Gloves** (when handling chemicals)
            - **Closed-toe shoes**

            Avoid:

            - Loose clothing
            - Dangling jewelry
            - Open hair (tie it back!)

            ✅ **Tip:** Dress for safety, not fashion, in the lab!
            """)

            # Section 3
            st.subheader("☣️ 3. Handle Chemicals with Care")
            st.markdown("""
            - Always read labels and **Safety Data Sheets (SDS)** before using a chemical.
            - Never taste or directly inhale chemicals.
            - Use fume hoods for volatile or hazardous substances.
            - Label all containers clearly.

            ✅ **Tip:** When diluting acids, always **add acid to water** — never the reverse.
            """)

            # Section 4
            st.subheader("⚗️ 4. Use Equipment Correctly")
            st.markdown("""
            - Inspect glassware for cracks before using.
            - Never use damaged or malfunctioning equipment.
            - Turn off burners and electrical devices after use.
            - Do not leave experiments unattended.

            ✅ **Tip:** Always follow the instructor’s or manual’s steps precisely.
            """)

            # Section 5
            st.subheader("🧼 5. Maintain Cleanliness")
            st.markdown("""
            - Clean spills immediately with proper materials.
            - Dispose of waste in designated containers — **not in the sink**.
            - Return tools and materials to their proper place after use.

            ✅ **Tip:** Leave the lab cleaner than you found it.
            """)

            # Section 6
            st.subheader("🚫 6. Know What NOT to Do")
            st.markdown("""
            - Never **eat, drink, or chew gum** in the lab.
            - Never mix chemicals for fun or without instruction.
            - Never touch your face, eyes, or mouth while working.

            ✅ **Tip:** Treat everything in the lab as potentially hazardous.
            """)

            # Section 7
            st.subheader("🧘 7. Stay Focused and Calm")
            st.markdown("""
            - Avoid horseplay or distractions.
            - If an accident occurs, **stay calm** and notify the supervisor immediately.
            - Report all injuries — even minor ones.

            ✅ **Tip:** Always prioritize safety over speed.
            """)

            # Final thoughts
            st.markdown("---")
            st.header("🌟 Final Thoughts")
            st.markdown("""
            Laboratory safety is **everyone’s responsibility**. By following these guidelines, you protect yourself, your peers, and the integrity of your experiments.

            Always stay alert, be respectful of lab rules, and treat chemicals and equipment with care. Remember — the best scientists are not only curious, but cautious.

            🧠 *“Safety isn’t expensive, it’s priceless.”*
            """)

        elif selected == "Career Guides":
            # Page title
            st.title("🧬 Career Guide: What Can You Do With a Degree in Chemistry?")

            st.markdown("""
            Chemistry isn’t just a subject in textbooks — it’s the foundation of numerous exciting and impactful careers.

            Whether you love working in a lab, want to explore environmental science, or dream of inventing new materials or medicines, **chemistry opens doors to diverse fields**. Here's a guide to help you discover your path!
            """)

            # Section: Why Choose Chemistry?
            st.subheader("🌟 Why Choose a Career in Chemistry?")
            st.markdown("""
            - It’s a **central science** connecting biology, physics, environmental science, and engineering.
            - Chemists contribute to **medicine, technology, sustainability, food safety**, and much more.
            - You can work in **academia, industry, government**, or **entrepreneurship**.

            Chemistry careers are **well-paying, intellectually rewarding**, and often contribute to real-world solutions.
            """)

            # Section: Top Career Options
            st.subheader("🧪 Top Career Options in Chemistry")

            st.markdown("""
            ### 1. **Analytical Chemist**
            - Role: Analyze substances to determine their composition.
            - Where: Pharma companies, food labs, forensic labs.
            - Skills: Spectroscopy, chromatography, critical thinking.

            ### 2. **Research Scientist**
            - Role: Discover new materials, reactions, or methods.
            - Where: Universities, R&D centers, biotech firms.
            - Skills: Experimental design, scientific writing, data analysis.

            ### 3. **Pharmaceutical Chemist**
            - Role: Develop and test new drugs.
            - Where: Biotech and pharmaceutical industries.
            - Skills: Organic chemistry, pharmacology, lab techniques.

            ### 4. **Forensic Scientist**
            - Role: Analyze crime scene evidence.
            - Where: Police labs, crime investigation units.
            - Skills: Toxicology, DNA analysis, analytical tools.

            ### 5. **Environmental Chemist**
            - Role: Monitor pollutants, develop green solutions.
            - Where: Government, environmental NGOs, consultancies.
            - Skills: Ecology, water/soil testing, sustainability knowledge.

            ### 6. **Chemical Engineer**
            - Role: Design chemical processes for manufacturing.
            - Where: Oil, gas, plastics, and food industries.
            - Skills: Process design, thermodynamics, simulation tools.

            ### 7. **Teaching & Academia**
            - Role: Teach or mentor the next generation of scientists.
            - Where: Schools, colleges, universities.
            - Skills: Communication, subject expertise, patience.

            ### 8. **Cosmetic Chemist**
            - Role: Formulate skincare, beauty, and personal care products.
            - Where: Cosmetic companies, R&D labs.
            - Skills: Formulation chemistry, microbiology, product testing.

            ### 9. **Toxicologist**
            - Role: Study the effects of chemicals on humans and the environment.
            - Where: Regulatory bodies, pharma, agriculture.
            - Skills: Risk assessment, biology, safety analysis.

            ### 10. **Science Communicator**
            - Role: Make science engaging for the public.
            - Where: Media, education platforms, museums.
            - Skills: Writing, video creation, simplifying complex topics.
            """)

            # Section: Skills You Need
            st.subheader("🧠 Essential Skills for Chemistry Careers")
            st.markdown("""
            - Laboratory techniques
            - Data analysis and software (Excel, Python, ChemDraw, etc.)
            - Research & scientific writing
            - Attention to detail
            - Problem-solving and teamwork
            - Ethical practices and safety

            ✅ **Tip:** Gain hands-on experience through internships or final-year projects.
            """)

            # Section: Educational Pathways
            st.subheader("🎓 What to Study?")
            st.markdown("""
            - **Bachelor's in Chemistry**: Foundation for lab roles and industry jobs.
            - **Master’s or PhD**: Required for advanced research, academia, or niche fields.
            - **Certifications**: OSHA lab safety, data analysis, or regulatory courses are helpful extras.

            📘 **Related Fields:** Biochemistry, Chemical Engineering, Environmental Science, Pharmacy.
            """)

            # Section: Final Advice
            st.markdown("---")
            st.header("💬 Final Advice")
            st.markdown("""
            - Be curious, not just about chemistry — but about where it fits in the world.
            - Network with professionals, attend science fairs, and stay updated on research.
            - Combine your chemistry background with other passions — like AI, law, business, or medicine.

            A degree in chemistry doesn’t limit you — **it empowers you** to pursue hundreds of directions. 🔬🌍

            💡 *“The country which is ahead in chemistry will be ahead in wealth and power.”* – William Ramsay
            """)

if __name__ == "__main__":
    display_blogs()