# Import packages needed
from util.load_packages import st, plt, np



# Load data
path = 'C:/Users/josef/Documents/GitHub/Master_CAU/Semester_3/Risk Management/Risk_App/data/DAX_index.csv'
data = np.genfromtxt(path, usecols=(1), delimiter=",", skip_header=1)

# Page configuration
st.set_page_config(page_title="Risk Management")

# Sidebar Navigation
st.sidebar.header("Navigate to Sections")
section = st.sidebar.radio("Select a section:", 
                           ["Introduction", "Methods"])




# Main page content
if section == "Introduction":
    st.title("Introduction to Risk Management")
    st.write("""
    Risk management is the process of identifying, assessing, and controlling threats to an organization's capital and earnings.
    It involves a systematic approach to managing risks to achieve the objectives of the organization.
    """)

    Subpage navigation within the Introduction section
    intro_section = st.radio("Select subpage:", ["Overview", "Risk Types", "Risk Management Process"])

    if intro_section == "Overview":
        st.write("The overview section explains the fundamental concepts of risk management.")
    elif intro_section == "Risk Types":
        st.write("This section describes the different types of risks, such as financial, operational, and strategic risks.")
    elif intro_section == "Risk Management Process":
        st.write("The risk management process includes risk identification, assessment, control, and monitoring.")









elif section == "Methods":
    st.title("Risk Management Methods")

    # Subpage navigation for Methods
    method_section = st.radio("Choose a method:", ["Normal Distribution", "Poisson Distribution", "Other Methods"])

    if method_section == "Normal Distribution":
        st.write("""
        The normal distribution is a fundamental concept in statistics and is widely used in risk management.
        It is often used to represent real-valued random variables whose distributions are not known.
        """)
        # Example of using the loaded data in a plot
        fig, ax = plt.subplots()
        ax.hist(data, bins=50, density=True)
        ax.set_title("Histogram of DAX Index")
        st.pyplot(fig)

    elif method_section == "Poisson Distribution":
        st.write("""
        The Poisson distribution models the number of events in fixed intervals of time or space.
        It's often used to model rare events or events with low probability.
        """)
        # Example of generating and plotting Poisson distribution
        poisson_data = np.random.poisson(lam=5, size=1000)
        fig, ax = plt.subplots()
        ax.hist(poisson_data, bins=30, density=True)
        ax.set_title("Poisson Distribution (λ=5)")
        st.pyplot(fig)

    elif method_section == "Other Methods":
        st.write("""
        This section discusses other methods used in risk management, such as Monte Carlo simulations or Value at Risk (VaR).
        """)
    
    
    
    #####
    
    
    
    
    
    
    
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
