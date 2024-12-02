## Import packages --
from util.load_packages import st



st.title("Risk Management Methods")


import pages.risk_measure.standard_deviation




st.write("""
Here are the risk management methods you can explore:
- Normal Distribution
- Poisson Distribution
- Other Methods
""")


# Sidebar Navigation
st.sidebar.header("Navigate to Sections")
section = st.sidebar.radio("Select a section:", 
                           ["Standard Deviation",
                            "Value at Risk",
                            "Expected Shortfall"])

# Main page content
if section == "Standard Deviation":
    st.title("Introduction to Risk Management")
    st.write("""
    Risk management is the process of identifying, assessing, and controlling threats to an organization's capital and earnings.
    It involves a systematic approach to managing risks to achieve the objectives of the organization.
    """)

    # Subpage navigation within the Introduction section
    import pages.risk_measure.standard_deviation




method_section = st.radio("Choose a method:", ["Normal Distribution", "Poisson Distribution", "Other Methods"])

# Dynamically import and execute the content based on selection
if method_section == "Standard Deviation":
    # Dynamically load normal_distribution.py content
    import pages.risk_measure.standard_deviation
elif method_section == "VaR":
    # Dynamically load poisson_distribution.py content
    import pages.risk_measure.VaR
elif method_section == "ES":
    # Dynamically load other_methods.py content
    import pages.risk_measure.ES
