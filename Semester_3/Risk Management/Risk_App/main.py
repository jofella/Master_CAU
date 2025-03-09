from util.load_packages import st  # Assuming `st` is imported from your utility script

# --- Head ---
st.set_page_config(
    page_title="Multipage Application",  # Page title for the browser tab
)

# --- Body ---
st.title("Main Page")  # Title for the page

# Main content
st.write("""
         - The general idea of this dashboard is to provide an interactive learning ressource. It combines both, 
         my theoretical/practical understanding about Risk Management and skillset in streamlit.
Expenses as material breeding insisted building to in. Continual so distrusts pronounce by 
unwilling listening. Thing do taste on we manor. Him had wound use found hoped. Of distrusts 
immediate enjoyment curiosity do. Marianne numerous saw thoughts the humoured.

And sir dare view but over man. So at within mr to simple assure. Mr disposing continued it 
offending arranging in we. Extremity as if breakfast agreement. Off now mistress provided out
horrible opinions. Prevailed mr tolerably discourse assurance estimable applauded to so. Him
everything melancholy uncommonly but solicitude inhabiting projection off. Connection stimulated
estimating excellence an to impression.
""")

# Add a separator line
st.write("---")

# File uploader for users to upload data
st.write("Upload your data here:")
uploaded_file = st.file_uploader("Choose a CSV file", type="csv")

if uploaded_file is not None:
    # If a file is uploaded, show its contents
    import pandas as pd
    df = pd.read_csv(uploaded_file)
    st.write(df.head())  # Display the first few rows of the uploaded data
