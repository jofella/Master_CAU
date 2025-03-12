from util.load_packages import st, np, pd, os, px, go, stats

st.title("📉 Losses")
st.markdown("""
Why do we need risk management? Financial institutions like banks are subject to losses. Extreme
large losses may lead to bancruptcy and may also threaten third parties. In order to prevent this
those institutions accumulate <b>buffer capital</b>. Three main questions arise: <br>

- How to quantify risk? <br>
- How to measure risk? <br>
- What capital reserve is needed in view of this risk? <br>

Following passages will introduce the basics for capturing <b>market risk</b>.
""", unsafe_allow_html=True)


st.write("---")

# ** Load DAX Companies Data **
@st.cache_data
def load_dax_data(path):
    """Load DAX stock data from CSV."""
    return np.genfromtxt(path, usecols=(1,2,3,4,5), delimiter=",", skip_header=1)

path = r'C:\Users\josef\Documents\GitHub\Master_CAU\Semester_3\Risk Management\Risk_App\data\DAX_companies.csv'
data_dax_comp = load_dax_data(path)



# ** Loss operator **
st.header("1. Loss Operator")

st.markdown("""
### 1.1. Understanding Risk Factors and the Loss Operator  

When looking at a (stock) portfolio in terms of risk, we chose **log stock prices** as  
<b>risk factors</b>. The risk factors are defined as:  

$$Z_{n,i} := \log(S_{n,i})$$  

Why risk factors? The idea is to look at the **key drivers of uncertainty** in the financial portfolio.  

Based on these, the **risk factor changes** are computed as:  

$$X_{n+1} = Z_{n+1} - Z_n$$  

These values serve as inputs for our **loss operator**, a function that calculates portfolio losses.  
The **advantage** of this approach is that we **separate risk factors from the portfolio structure**,  
making the modeling process **more flexible and less complicated**.  

---

### 1.2. The Loss Operator  

The loss in period \( n+1 \) is a function of \( X_{n+1} \) and known quantities at \( t_n \).  
Specifically, we define the **loss operator** as:

$$
L_{n+1} = -(V_{n+1} - V_n) = - f_{n+1}(Z_n + X_{n+1}) + f_n(Z_n) =: \ell_{[n]}(X_{n+1})
$$

The function L is called the **loss operator**, which **randomly changes** over time. <br>

""", unsafe_allow_html=True)

st.write("---")

st.write("""
Below you can see set your desired portfolio weights and see how it affects both: portfolio value and losses.
""")

col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    w1 = st.number_input("Asset 1", min_value=0, max_value=30, value=4, step=1)
with col2:
    w2 = st.number_input("Asset 2", min_value=0, max_value=30, value=8, step=1)
with col3:
    w3 = st.number_input("Asset 3", min_value=0, max_value=30, value=15, step=1)
with col4:
    w4 = st.number_input("Asset 4", min_value=0, max_value=30, value=16, step=1)
with col5:
    w5 = st.number_input("Asset 5", min_value=0, max_value=30, value=23, step=1)

alpha_weights = np.array([w1, w2, w3, w4, w5])


# ** Computation **

# Risk factors
@st.cache_data
def compute_risk_factors(data):
    """Compute log stock prices (Z_n) and log-returns (X_n)."""
    Z_n = np.log(data)
    X_n = np.diff(Z_n, axis=0)  # Compute log-return changes
    return Z_n, X_n

Z_n, X_n = compute_risk_factors(data_dax_comp)

# Weights
weighted_port = alpha_weights * data_dax_comp

# Compute Portfolio Value (V_n)
V_n = np.dot(np.exp(Z_n), alpha_weights)

# Losses
@st.cache_data
def compute_nonlinear_losses(X_n, alpha_weights, data):
    """Compute the nonlinear portfolio losses."""
    weighted_port = alpha_weights * data
    
    def l(n, x):
        return -np.dot(weighted_port[n, :], np.exp(x[n, :]) - 1)
    
    return np.array([l(n, X_n) for n in range(len(X_n))])

@st.cache_data
def compute_linearized_losses(X_n, alpha_weights, data):
    """Compute the linearized (first-order Taylor) portfolio losses."""
    weighted_port = alpha_weights * data
    
    def l_delta(n, x):
        return -np.dot(weighted_port[n, :], x[n, :])
    
    return np.array([l_delta(n, X_n) for n in range(len(X_n))])

losses = compute_nonlinear_losses(X_n, alpha_weights, data_dax_comp)
delta_losses = compute_linearized_losses(X_n, alpha_weights, data_dax_comp)

# Convert to DataFrames
loss_df = pd.DataFrame({"Time": np.arange(len(losses)), "Losses": losses})
delta_loss_df = pd.DataFrame({"Time": np.arange(len(delta_losses)), "Losses": delta_losses})




# ** DAX portfolio **
fig = go.Figure()

fig.add_trace(go.Scatter(
    x=np.arange(len(V_n)),
    y=V_n,
    mode="lines",
    name="Simulated DAX Portfolio",
    line=dict(color="light blue")
))

fig.update_layout(
    title="DAX 5-stock Portfolio from 2000-Today",
    xaxis_title="Time (Days)",
    yaxis_title="Portfolio Value",
    xaxis=dict(showgrid=True),
    yaxis=dict(showgrid=True)
)

st.plotly_chart(fig)



# ** Losses **
fig = px.line(loss_df, x="Time", y="Losses",
              title="Portfolio Losses Over Time",
              labels={"Losses": "Loss", "Time": "Time (Days)"},
              line_shape="linear")

st.plotly_chart(fig)


st.write("---")

st.markdown("""
### 1.3. Linearized Loss Operator  

One major challenge with the **loss operator** is that it is generally **nonlinear** due to the function \( f_n \).  
Nonlinearities can make risk modeling more complex, requiring advanced numerical methods for estimation.  

To simplify computations and make risk estimation more tractable, we apply a **first-order Taylor expansion** to approximate losses **linearly**.  
This approximation allows us to efficiently analyze **small risk factor changes** while maintaining accuracy.

""", unsafe_allow_html=True)



@st.cache_data
def plot_losses_chart(loss_df, delta_loss_df):
    """Generates a Plotly figure comparing nonlinear & linearized losses."""
    fig = go.Figure()

    # Nonlinear Losses
    fig.add_trace(go.Scatter(
        x=loss_df["Time"], 
        y=loss_df["Losses"], 
        mode="lines", 
        name="Nonlinear Losses", 
        line=dict(color="lightblue")
    ))

    # Linearized Losses
    fig.add_trace(go.Scatter(
        x=delta_loss_df["Time"], 
        y=delta_loss_df["Losses"], 
        mode="lines", 
        name="Linearized Losses", 
        line=dict(color="crimson", dash="dash")
    ))

    fig.update_layout(
        title="Comparison of Nonlinear and Linearized Portfolio Losses",
        xaxis_title="Time (Days)",
        yaxis_title="Loss",
        legend=dict(x=0, y=1),
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True)
    )

    return fig

fig = plot_losses_chart(loss_df, delta_loss_df)
st.plotly_chart(fig)


st.write("---")

