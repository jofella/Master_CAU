# All commonly used function like load, etc.
from util.load_packages import pd, np, st


# 1. Get log-returns
def get_log_returns(data):
    return np.diff(np.log(data))


# 2. Load (single-stock) data
@st.cache_data
def load_single_stock_data(path):
    """Load single-stock data from CSV."""
    return np.genfromtxt(path, usecols=(1), delimiter=",", skip_header=1)

