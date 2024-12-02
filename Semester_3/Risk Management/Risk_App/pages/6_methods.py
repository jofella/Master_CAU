import streamlit as st

st.title("Risk Management Methods")

st.write("""
Here are the risk management methods you can explore:
- Normal Distribution
- Poisson Distribution
- Other Methods
""")


# Sidebar Navigation
st.sidebar.header("Navigate to Sections")
section = st.sidebar.radio("Select a section:", 
                           ["Normal", "Poisson"])
# Main page content
if section == "Normal":
    st.title("Introduction to Risk Management")
    st.write("""
    Risk management is the process of identifying, assessing, and controlling threats to an organization's capital and earnings.
    It involves a systematic approach to managing risks to achieve the objectives of the organization.
    """)

    # Subpage navigation within the Introduction section
    import pages.methods.normal_distribution



method_section = st.radio("Choose a method:", ["Normal Distribution", "Poisson Distribution", "Other Methods"])

# Dynamically import and execute the content based on selection
if method_section == "Normal Distribution":
    # Dynamically load normal_distribution.py content
    import pages.methods.normal_distribution
elif method_section == "Poisson Distribution":
    # Dynamically load poisson_distribution.py content
    import pages.methods.poisson_distribution
elif method_section == "Other Methods":
    # Dynamically load other_methods.py content
    import pages.methods.other_methods
