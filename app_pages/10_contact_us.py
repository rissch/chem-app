import streamlit as st
import smtplib
from email.message import EmailMessage
import os

st.title("📩 Contact Us")

with st.form("contact_form"):
    name = st.text_input("Name *", max_chars=50)
    email = st.text_input("Email *", max_chars=100)
    message = st.text_area("Message", height=150)
    submitted = st.form_submit_button("Send")

if submitted:
    if not name.strip():
        st.error("Please enter your name.")
    elif not email.strip():
        st.error("Please enter your email.")
    else:
        # Prepare email
        SMTP_USERNAME = os.getenv("social.personalcalc@gmail.com")           # Your Gmail address
        SMTP_PASSWORD = os.getenv("O!x%ce#E^s9eLHLwkPgWtG8A")    # Your Gmail app password
        email_message = EmailMessage()
        email_message['Subject'] = f"New Contact Form Submission from {name}"
        email_message['From'] = SMTP_USERNAME             # Sender email must be your verified SMTP email
        email_message['To'] = "hello@molaritycalculator.net"
        email_message['Reply-To'] = email                  # User's email for reply

        body = f"Name: {name}\nEmail: {email}\n\nMessage:\n{message if message.strip() else '(No message provided)'}"
        email_message.set_content(body)

        # SMTP server configuration
        SMTP_SERVER = "smtp.gmail.com"
        SMTP_PORT = 587

        try:
            with smtplib.SMTP(SMTP_SERVER, SMTP_PORT) as smtp:
                smtp.starttls()
                smtp.login(SMTP_USERNAME, SMTP_PASSWORD)
                smtp.send_message(email_message)
            st.success("Thank you! Your message has been sent.")
        except Exception as e:
            st.error(f"Oops! Something went wrong while sending your message: {e}")
