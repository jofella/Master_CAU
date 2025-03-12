from util.load_packages import st, np, pd, os, px, go, stats
from util.data_utils import get_log_returns

# Normal Distribution section content
st.markdown("""
### 1. Standard Deviation  

A **traditional approach** to measuring risk is based on the **standard deviation of the loss**.  
It quantifies the typical deviation of losses from their expected value, providing a **measure of risk volatility**:

$$
\\varrho = c \\sqrt{\\operatorname{Var}_n(L_{n+1})} = c \\sqrt{ \\int x^2 P^L(dx) - \\left( \\int x P^L(dx) \\right)^2 }
$$

where \( c > 0 \) is a constant factor, possibly adjusted for the mean by adding \( E_n(L_{n+1}) \) if necessary.


### Benefits of Standard Deviation as a Risk Measure
✔ **Simple & Easy to Estimate** – Standard deviation is a widely used and well-understood metric.  
✔ **Symmetric Risk Interpretation** – It treats **profits and losses equally**, making it useful for normal distributions.  
✔ **Provides a General Risk Indicator** – Helps assess **volatility and overall dispersion of losses**.  


### Limitations of Standard Deviation for Risk Assessment  
🚨 **Fails to Capture Tail Risk** – It does not focus on extreme losses, which are **critical for risk management**.  
🚨 **Assumes a Symmetric Distribution** – In financial markets, loss distributions are often **skewed** and **fat-tailed**.  
🚨 **Not Ideal for Heavy-Tailed Risks** – If losses have **infinite variance**, standard deviation **is not well-defined**.  

### Key Takeaway  
While **standard deviation is useful as a general risk measure**, it may be **insufficient** for capturing extreme financial risks.  
For more robust risk assessment, **Value at Risk (VaR) and Expected Shortfall (ES)** are often preferred.

""", unsafe_allow_html=True)

st.write("---")


# **1.Load data **
data = st.session_state.data


# ** 2.Computation of sd**
def rho(n, c, mu, sigma):
    cond_mean = n * (1 - np.exp(mu+ sigma**2 / 2))
    cond_var = n**2 * (np.exp(sigma**2) - 1) * np.exp(2*mu + sigma**2)
    return cond_mean + c * np.sqrt(cond_var)
c = 1.64


st.write(f"Pre-defined risk-appetite is: {c}")

# lr
returns_dax = get_log_returns(data)

# Losses
losses = -np.diff(data)

# Get parameters
mu = np.mean(returns_dax)
sigma = np.std(returns_dax)

# compute standard deviation
sd = np.empty(len(data))

for i in range(len(data)):
    sd[i] = rho(data[i], c, mu, sigma)


# Transform to df
loss_df = pd.DataFrame({"Index": range(len(losses)), "Losses": losses})
sd_df = pd.DataFrame({"Index": range(len(data)), "sd": sd})



# ** 2.Visualization **
@st.cache_data
def plot_sd(loss_df, sd_df):
    """Generates a Plotly figure comparing losses and sd."""
    fig = go.Figure()

    # Losses
    fig.add_trace(go.Scatter(
        x=loss_df.index, 
        y=loss_df["Losses"], 
        mode="lines", 
        name="Losses", 
        line=dict(color="lightblue")
    ))

    # sd
    fig.add_trace(go.Scatter(
        x=sd_df.index, 
        y=sd_df["sd"], 
        mode="lines", 
        name="Standard Deviation", 
        line=dict(color="crimson")
    ))

    fig.update_layout(
        # title="Standard Deviation",
        xaxis_title="Time (Days)",
        yaxis_title="Loss",
        legend=dict(x=0, y=1),
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True)
    )

    return fig

fig = plot_sd(loss_df, sd_df)
st.plotly_chart(fig)


st.write("---")