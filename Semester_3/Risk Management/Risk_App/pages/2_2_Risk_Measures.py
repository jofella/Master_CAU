from util.load_packages import st, np, pd, os, px, go, stats


st.title("📏 Measuring Risk")

st.markdown("""
We have now established a method to quantify our risk: **Losses**.  
The next step is to measure risk and ultimately determine the **required buffer capital**.  
""", unsafe_allow_html=True)


# ** 1.Load data (once and reference on it) **
if "data" not in st.session_state:
    st.session_state.data = pd.read_csv(r'C:\Users\josef\Documents\GitHub\Master_CAU\Semester_3\Risk Management\Risk_App\data\DAX_companies.csv')  # Replace with actual path



method_section = st.radio("Choose a risk measure:", ["Standard Deviation", "Value at Risk", "Expected Shortfall"])

# Dynamically import and execute the content based on selection
if method_section == "Standard Deviation":
    # Dynamically load normal_distribution.py content
    import pages.risk_measure.standard_deviation
elif method_section == "Value at Risk":
    # Dynamically load poisson_distribution.py content
    import pages.risk_measure.VaR
elif method_section == "Expected Shortfall":
    # Dynamically load other_methods.py content
    import pages.risk_measure.ES
