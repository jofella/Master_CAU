## Import packages --
from util.load_packages import st, np, pd, os, px, go, stats
from util.load_packages import st



### Body ---
st.title("Risk Management Methods")
import pages.risk_measure.standard_deviation
import pages.risk_measure.VaR
import pages.risk_measure.ES

st.header("1. Loss operator")

st.markdown("""Why does it make sense to compute the losses from the risk factor changes instead of the time series directly?
            By splitting the calculation of the losses into the risk factor changes and the loss operator we get the advantage that we can now model both separately. The loss operator depends on the structure of the portfolio and the chosen risk factors. If you now want to incorporate a new risk factor or the structure of your portfolio changes then you only have to adjust the loss operator slightly. If you want to take a new model for your risk factor changes, then you only have to change these. For example the risk factor changes in this exercise are the log-returns. 
            If you want to do a Monte-Carlo simulation you can now choose any fitting distribution without any thought about the loss operator.
            """)


