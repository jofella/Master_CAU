import numpy as np
import plotly.graph_objects as go
import streamlit as st


st.title("📏 Measuring Risk")
st.markdown("""
Introduction text ...
""")


st.write("---")


# ** Load DAX Companies Data **
path = r'C:\Users\josef\Documents\GitHub\Master_CAU\Semester_3\Risk Management\Risk_App\data\DAX_companies.csv'
data_dax_comp = np.genfromtxt(path, usecols=(1,2,3,4,5), delimiter=",", skip_header=1)

# ** Compute Risk Factors (Z_n) **
Z_n = np.log(data_dax_comp)


# ** Section: Simulate Portfolio Process **
st.header("Simulated DAX Portfolio")
st.markdown("""
This section visualizes the **DAX portfolio performance** over time. 
The portfolio is constructed using **log returns** and assigned weights to individual stocks.
""")


col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    w1 = st.slider("Asset 1", min_value=0, max_value=30, value=4, step=1)
with col2:
    w2 = st.slider("Asset 2", min_value=0, max_value=30, value=8, step=1)
with col3:
    w3 = st.slider("Asset 3", min_value=0, max_value=30, value=15, step=1)
with col4:
    w4 = st.slider("Asset 4", min_value=0, max_value=30, value=16, step=1)
with col5:
    w5 = st.slider("Asset 5", min_value=0, max_value=30, value=23, step=1)

alpha_weights = np.array([w1, w2, w3, w4, w5])



# ** Compute Portfolio Value (V_n) **
V_n = np.dot(np.exp(Z_n), alpha_weights)



# ** Plot with Plotly **
fig = go.Figure()

fig.add_trace(go.Scatter(
    x=np.arange(len(V_n)),
    y=V_n,
    mode="lines",
    name="Simulated DAX Portfolio",
    line=dict(color="light blue")
))

fig.update_layout(
    title="DAX Portfolio from 2000-Today",
    xaxis_title="Time (Days)",
    yaxis_title="Portfolio Value",
    xaxis=dict(showgrid=True),  # Enable x-axis grid
    yaxis=dict(showgrid=True)   # Enable y-axis grid
)

# ** Display in Streamlit **
st.plotly_chart(fig)

st.markdown("""
### 📝 Interpretation:
- The **simulated portfolio** follows an **exponential growth pattern**, consistent with a compounding return process.
- Since we apply **log transformations and exponentiation**, the result reflects a **multiplicative** return process.
- The weights `[4, 8, 15, 16, 23]` contribute differently to portfolio fluctuations.
""")



# Calculate Losses

# Set parameters
alpha_weights = [4, 8, 15, 16, 23]
weighted_port = alpha_weights * data_dax_comp

# Function loss operator
def l(n, x):
    return -np.dot(weighted_port[n,:], np.exp(x[n,:])-1)



# Get Losses (with function)
baba_matrix = np.zeros(len(X_n_1))

for i in range(len(X_n_1)):
    baba_matrix[i] = l(i, X_n_1)
    
baba_matrix


# Delta loss operator


