from util.load_packages import st, np, pd, os, px, go, stats
from util.data_utils import get_log_returns


st.markdown("""
### 2.2. Value at Risk


""", unsafe_allow_html=True)


# **1.Load data **
data = st.session_state.data


# **2.Computation **
returns_dax = get_log_returns(data)

losses = -np.diff(data)
mu = np.mean(returns_dax)
sigma = np.std(returns_dax)



# **3.VaR vs Lin. VaR **



st.write("""
- first the function
- secondly backtesting the VaR
""")

st.write("---")


st.write("---")