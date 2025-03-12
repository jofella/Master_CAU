from util.load_packages import st, np, pd, os, px, go, stats
from util.data_utils import load_single_stock_data

st.title("📏 Measuring Risk")

st.markdown("""
We have now established a method to quantify our risk: **Losses**.  
The next step is to measure risk and ultimately determine the **required buffer capital**.  
""", unsafe_allow_html=True)

# ** 1.Load data (once and reference on it) **
path = r'C:\Users\josef\Documents\GitHub\Master_CAU\Semester_3\Risk Management\Risk_App\data\DAX_companies.csv'

# Load data and store it in session state
if "data_dax_comp" not in st.session_state:
    st.session_state.data = load_single_stock_data(path)



# ** 2.Select desired risk measure **
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
