from util.load_packages import st, np, pd, os, px, go, stats

st.title("📏 Measuring Risk")
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
path = r'C:\Users\josef\Documents\GitHub\Master_CAU\Semester_3\Risk Management\Risk_App\data\DAX_companies.csv'
data_dax_comp = np.genfromtxt(path, usecols=(1,2,3,4,5), delimiter=",", skip_header=1)


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
# Compute Risk Factors (Z_n) 
Z_n = np.log(data_dax_comp)

# Compute Risk Factor Changes (X_n)
X_n = np.diff(Z_n, axis=0)  # Compute log-return changes

weighted_port = alpha_weights * data_dax_comp

# Function loss operator
def l(n, x):
    return -np.dot(weighted_port[n,:], np.exp(x[n,:])-1)

# Compute Losses Over Time
losses = np.array([l(n, X_n) for n in range(len(X_n))])

# Convert to DataFrame
V_n = np.dot(np.exp(Z_n), alpha_weights)
loss_df = pd.DataFrame({"Time": np.arange(len(losses)), "Losses": losses})







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

# ** Display in Streamlit **
st.plotly_chart(fig)




st.write("---")


st.markdown("""
### 1.3. Linearlized Loss operator

hi



""", unsafe_allow_html=True)
