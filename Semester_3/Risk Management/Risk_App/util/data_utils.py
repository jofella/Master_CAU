# All commonly used function like load, etc.
from util.load_packages import pd, np, st


# 1. Get log-returns
def get_log_returns(data):
    return np.diff(np.log(data))

# 2. Loss operator