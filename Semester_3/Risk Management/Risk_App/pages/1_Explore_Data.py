from util.load_packages import st, np, pd, os, px, go, stats, plt
from util.data_utils import get_log_returns

# Page content for Explore Data
st.title("Explore Data")

st.write("""
This page allows you to upload and explore your data.
""")

st.write("The first step is to get the returns of ts and look at the data visually:")

## Load data ---
st.write(f"Looking for file at: {os.path.abspath('..Risk_App/data/DAX_index.csv')}")

path = r'C:\Users\josef\Documents\GitHub\Master_CAU\Semester_3\Risk Management\Risk_App\data\DAX_index.csv'
data = np.genfromtxt(path, usecols=(1), delimiter=",", skip_header=1)



## Get log-returns --
lr = get_log_returns(data)

# Compute mean and standard deviation
mu = np.mean(lr)
sigma = np.std(lr)

# Generate normal distribution for overlay
x_rand_range = np.linspace(min(lr), max(lr), 100)
norm_pdf = stats.norm.pdf(x_rand_range, loc=mu, scale=sigma)

# Scale the normal PDF to match the histogram height
bin_width = (max(lr) - min(lr)) / 1000  # Adjust bins to match histogram
norm_pdf_scaled = norm_pdf * len(lr) * bin_width

# Convert empirical log returns to DataFrame
lr_df = pd.DataFrame({"Index": range(len(lr)), "Log Returns": lr})

# MC-simulation for synthetic log returns
simulated_r = np.random.normal(mu, sigma, len(lr))

# Convert simulated returns to DataFrame
sim_df = pd.DataFrame({"Index": range(len(simulated_r)), "Simulated Log Returns": simulated_r})



## Line Chart (Empirical vs. Simulated Log Returns) --
fig = go.Figure()

# Line for actual log returns
fig.add_trace(go.Scatter(
    x=lr_df["Index"], 
    y=lr_df["Log Returns"], 
    mode="lines", 
    name="Empirical Log Returns",
    line=dict(color="light blue")
))

# Line for simulated log returns
fig.add_trace(go.Scatter(
    x=sim_df["Index"], 
    y=sim_df["Simulated Log Returns"], 
    mode="lines", 
    name="Simulated Log Returns",
    line=dict(color="red", dash="dash")
))

fig.update_layout(
    title="Empirical vs Simulated Log Returns",
    xaxis_title="Index",
    yaxis_title="Log Returns",
    legend=dict(x=0, y=1)
)

st.plotly_chart(fig)


st.write("---")

## 📊 **Histogram with Normal Distribution Overlay**
fig = px.histogram(lr_df, x="Log Returns",
                   nbins=1000,
                   title="Plotly Histogram of Log Returns",
                   histnorm='probability density')

# Overlay Normal PDF
fig.add_trace(go.Scatter(
    x=x_rand_range,
    y=norm_pdf_scaled,  # Use scaled PDF
    mode="lines",
    name="Normal PDF",
    line=dict(color="red", width=2)
))

# Display in Streamlit
st.plotly_chart(fig)

st.write("---")
