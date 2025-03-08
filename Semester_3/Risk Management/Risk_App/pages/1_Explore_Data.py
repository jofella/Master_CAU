from util.load_packages import st, np, pd, plt, os, px
from util.data_utils import get_log_returns

# Page content for Explore Data
st.title("Explore Data")

st.write("""
This page allows you to upload and explore your data.
""")


## Load data ---
st.write(f"Looking for file at: {os.path.abspath('..Risk_App/data/DAX_index.csv')}")

path = path = r'C:\Users\josef\Documents\GitHub\Master_CAU\Semester_3\Risk Management\Risk_App\data\DAX_index.csv'
data = np.genfromtxt(path,
                     usecols=(1),
                     delimiter=",",
                     skip_header=1)


## Get log-returns --
lr = get_log_returns(data)


## Get plots --

# Line chart
df = pd.DataFrame({"Index": range(len(lr)), "Log Returns": lr})
fig = px.line(df, x="Index", y="Log Returns", title="Plotly Log Returns Line Chart")
st.plotly_chart(fig)


# Histogram
fig = px.histogram(df, x="Log Returns", nbins=30, title="Plotly Histogram of Log Returns")
st.plotly_chart(fig)