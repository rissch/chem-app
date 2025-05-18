import streamlit as st
from streamlit_option_menu import option_menu


def display_unit_converters():
    calc = option_menu(
        menu_title="Unit Converters",
        options=["Mass", "Volume", "Temperature", "Pressure", "Energy"],
        icons=["weight", "cup-straw", "thermometer", "speedometer", "battery-charging"],
        orientation="horizontal"
    )
    st.header(f"{calc} Converter")

    if calc == "Mass":
        st.subheader("⚖️ Mass Converter")
        mass = st.number_input("Enter mass", value=1.0)
        from_unit = st.selectbox("From", ["g", "kg", "mg"], index=1)
        to_unit = st.selectbox("To", ["g", "kg", "mg"], index=0)

        conversion = {"g": 1, "kg": 1000, "mg": 0.001}
        result = mass * conversion[from_unit] / conversion[to_unit]
        st.success(f"{mass} {from_unit} = {result:.4f} {to_unit}")

    elif calc == "Volume":
        st.subheader("🧪 Volume Converter")
        volume = st.number_input("Enter volume", value=1.0)
        from_unit = st.selectbox("From", ["L", "mL", "cm³"], index=0)
        to_unit = st.selectbox("To", ["L", "mL", "cm³"], index=1)

        conversion = {"L": 1000, "mL": 1, "cm³": 1}
        result = volume * conversion[from_unit] / conversion[to_unit]
        st.success(f"{volume} {from_unit} = {result:.4f} {to_unit}")

    elif calc == "Temperature":
        st.subheader("🌡️ Temperature Converter")
        temp = st.number_input("Enter temperature")
        from_unit = st.selectbox("From", ["Celsius", "Kelvin", "Fahrenheit"], index=0)
        to_unit = st.selectbox("To", ["Celsius", "Kelvin", "Fahrenheit"], index=2)

        def convert_temp(value, from_u, to_u):
            if from_u == to_u:
                return value
            if from_u == "Celsius":
                if to_u == "Kelvin": return value + 273.15
                elif to_u == "Fahrenheit": return (value * 9/5) + 32
            if from_u == "Kelvin":
                if to_u == "Celsius": return value - 273.15
                elif to_u == "Fahrenheit": return (value - 273.15) * 9/5 + 32
            if from_u == "Fahrenheit":
                if to_u == "Celsius": return (value - 32) * 5/9
                elif to_u == "Kelvin": return (value - 32) * 5/9 + 273.15

        result = convert_temp(temp, from_unit, to_unit)
        st.success(f"{temp} {from_unit} = {result:.2f} {to_unit}")

    elif calc == "Pressure":
        st.subheader("🧭 Pressure Converter")
        pressure = st.number_input("Enter pressure", value=1.0)
        from_unit = st.selectbox("From", ["atm", "Pa", "kPa", "mmHg", "bar"], index=1)
        to_unit = st.selectbox("To", ["atm", "Pa", "kPa", "mmHg", "bar"], index=2)

        conversion = {
            "atm": 101325,
            "Pa": 1,
            "kPa": 1000,
            "mmHg": 133.322,
            "bar": 100000
        }

        result = pressure * conversion[from_unit] / conversion[to_unit]
        st.success(f"{pressure} {from_unit} = {result:.3f} {to_unit}")

    elif calc == "Energy":
        st.subheader("⚡ Energy Converter")
        energy = st.number_input("Enter energy", value=1.0)
        from_unit = st.selectbox("From", ["J", "kJ", "cal", "kcal"], index=1)
        to_unit = st.selectbox("To", ["J", "kJ", "cal", "kcal"], index=3)

        conversion = {
            "J": 1,
            "kJ": 1000,
            "cal": 4.184,
            "kcal": 4184
        }

        result = energy * conversion[from_unit] / conversion[to_unit]
        st.success(f"{energy} {from_unit} = {result:.4f} {to_unit}")
