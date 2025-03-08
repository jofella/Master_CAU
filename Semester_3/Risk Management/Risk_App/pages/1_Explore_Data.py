from util.load_packages import st, np, pd, os, px, go, stats
from util.data_utils import get_log_returns


# ** 1.Introduction **
st.title("📊 Explore Data")

st.markdown("""
This dashboard helps analyze **log returns**, compare them with a **normal distribution**, 
and simulate returns using Monte Carlo methods.
""")

st.write("---")


# ** 2.Data Loading **
st.header("1. Load Data")
st.write("We will use historical data to compute log returns and compare them with simulations.")

st.write(f"Looking for file at: {os.path.abspath('..Risk_App/data/DAX_index.csv')}")

path = r'C:\Users\josef\Documents\GitHub\Master_CAU\Semester_3\Risk Management\Risk_App\data\DAX_index.csv'
data = np.genfromtxt(path, usecols=(1), delimiter=",", skip_header=1)


## Get log-returns --
lr = get_log_returns(data)

mu = np.mean(lr)
sigma = np.std(lr)

# Get norm. distr
x_rand_range = np.linspace(min(lr), max(lr), 100)
norm_pdf = stats.norm.pdf(x_rand_range, loc=mu, scale=sigma)

bin_width = (max(lr) - min(lr)) / 1000
norm_pdf_scaled = norm_pdf * len(lr) * bin_width # Scale normal PDF (match hist. height)

# MC-simulation for synthetic log returns
sim_returns = np.random.normal(mu, sigma, len(lr))

# Simulate price process
S_0 = data[0]
S_t = S_0 * np.exp(np.cumsum(sim_returns))

# Get DataFrames (plotly)
lr_df = pd.DataFrame({"Index": range(len(lr)), "Log Returns": lr})
sim_df = pd.DataFrame({"Index": range(len(sim_returns)), "Simulated Log Returns": sim_returns})
price_df = pd.DataFrame({"Index": range(len(data)), "Real DAX": data})
sim_price_df = pd.DataFrame({"Index": range(len(S_t)), "Simulated DAX": S_t})



st.write("---")


# ** 2.Key Metrics **
st.header("2. Summary Statistics")
st.write("Below the key moments of your log-retruns:")

mu = np.mean(lr)
sigma = np.std(lr)

col1, col2 = st.columns(2)
with col1:
    st.metric(label="Mean (μ)", value=f"{mu:.6f}")
with col2:
    st.metric(label="Standard Deviation (σ)", value=f"{sigma:.6f}")

st.markdown("These are used to simulate a normal pdf which is shown below. "
            "Here we use 'numpy', but any other statstical package can be used as well.")


st.write("---")


# ** 3.Visualizations **
st.header("3. Log Returns Over Time")

fig = go.Figure()

fig.add_trace(go.Scatter(
    x=lr_df["Index"], 
    y=lr_df["Log Returns"], 
    mode="lines", 
    name="Empirical Log Returns",
    line=dict(color="light blue")
))

# Overlay sim lr
fig.add_trace(go.Scatter(
    x=sim_df["Index"], 
    y=sim_df["Simulated Log Returns"], 
    mode="lines", 
    name="Simulated Log Returns",
    line=dict(color="red", dash="dash")
))

fig.update_layout(
    xaxis_title="Index",
    yaxis_title="Log Returns",
    legend=dict(x=0, y=1)
)

st.plotly_chart(fig)


st.markdown("""
    <p style="font-size:20px;">
    <b>Findings:</b><br>
    - Simulated returns appear to be more "ordered". Basically a stationary time series. <br>
    - The empirical returns are much more random in a sense, that we have large spikes not "covered" by the normal returns.<br>
    -> This indicates the good old story of the normal dist. <b> underestimating </b> extreme events. 
    </p>
    """, unsafe_allow_html=True)


st.write("---")


# **4.Histogram vs Normal Distribution**
st.header("4. Distribution of Log Returns")

num_bins = st.slider("Select Number of Bins:",
                     min_value=20,
                     max_value=500,
                     value=250,
                     step=10)

fig = px.histogram(lr_df, x="Log Returns",
                   nbins=num_bins,
                   histnorm='probability density'
                   )

# Overlay Normal PDF
fig.add_trace(go.Scatter(
    x=x_rand_range,
    y=norm_pdf_scaled,  # Use scaled PDF
    mode="lines",
    name="Normal PDF",
    line=dict(color="red", width=2)
))

st.plotly_chart(fig)

st.markdown("""
    <p style="font-size:20px;">
    <b>Findings:</b><br>
    - Similar story than in 3: Especially in the tails the normal distribution doesn't seem to be appropriate for modeling returns.
    The ones observed in the markets do behave differently.
    </p>
    """, unsafe_allow_html=True)



st.write("---")



# **5.Simulate Process**
st.header("5. Simulate Price Process")

fig = go.Figure()

# Plot Real DAX Prices
fig.add_trace(go.Scatter(
    x=price_df["Index"], 
    y=price_df["Real DAX"], 
    mode="lines",
    name="Real DAX",
    line=dict(color="light blue")
))

# Overlay Simulated Prices
fig.add_trace(go.Scatter(
    x=sim_price_df["Index"], 
    y=sim_price_df["Simulated DAX"], 
    mode="lines",
    name="Simulated DAX",
    line=dict(color="red", dash="dash")
))

fig.update_layout(
    xaxis_title="Index",
    yaxis_title="DAX Price",
    legend=dict(x=0, y=1)
)

st.plotly_chart(fig)


st.markdown("""
    <p style="font-size:20px;">
    <b>Findings:</b><br>
    - Attention: This will change with every re-load! <br>
    - Real DAX prices follow an upward drift. Its influenced by actual market conditions like eco. growth or inflation. <br>
    - Simulated DAX is much more volatile due to constant (stationary) volatilty. Much larger fluctuations, 
    beacuse here we dont account for different vola periods. <br>
    - Simulated process follows exponential growth path, massively overestimates vola 
    --> Lacks mean-reverting behavior. (Improvement: GARCH model)
    </p>
    """, unsafe_allow_html=True)


# st.write("---")

# # **5. Conclusions**
# st.header("5. Key Insights")

# st.markdown("""
# - The log returns appear **centered around zero** with some skewness.
# - The **simulated normal distribution** provides a baseline for risk comparison.
# - Future work: Try using **GARCH models** to estimate volatility! 🔥
# """)