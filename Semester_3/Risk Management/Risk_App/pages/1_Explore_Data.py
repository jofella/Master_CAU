from util.load_packages import st, np, plt, os

# Page content for Explore Data
st.title("Explore Data")

st.write("""
This page allows you to upload and explore your data.
""")

# See where its looking fro
print(f"Looking for file at: {os.path.abspath('../data/DAX_index.csv')}")


## Load data --
data = np.genfromtxt('.^data/DAX_index.csv',
                     usecols=(1), delimiter=",", skip_header=1) #./ --> in diesem Verzeichnis

fig, ax = plt.subplots()
ax.hist(data, bins=30, density=True)
st.pyplot(fig)

