import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Poisson Distribution section content
st.write("""
The Poisson Distribution models the number of events occurring within a fixed interval of time or space.
Here's some content specific to the Poisson Distribution.
""")

# Example: Plot a Poisson distribution
data = np.random.poisson(5, 1000)

fig, ax = plt.subplots()
ax.hist(data, bins=30, density=True)
st.pyplot(fig)
