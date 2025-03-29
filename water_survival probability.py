import MDAnalysis
from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

Model=3
Rep=0


# exponencial
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

# Analysis
universe = MDAnalysis.Universe("hCx46_"+str(Model)+"Ca.psf", "MD_500ns_"+str(Model)+"Ca_r"+str(Rep)+".dcd")
select = "byres name OH2 and cyzone 10.0 15.0 -15.0 protein"
sp = SP(universe, select, verbose=True)
sp.run(start=0, stop=499, tau_max=5)
tau_timeseries = sp.tau_timeseries
sp_timeseries = sp.sp_timeseries

# Adjust
popt, pcov = curve_fit(func, tau_timeseries, sp_timeseries)
trend_x = np.linspace(min(tau_timeseries), max(tau_timeseries), 100)
trend_y = func(trend_x, *popt)

# Graph
plt.figure(figsize=(16, 12))
plt.plot(tau_timeseries, sp_timeseries, 'bo', label='Original Data', markersize=10)
plt.plot(trend_x, trend_y, 'r-', label='Trendline: $f(x) = %5.3f e^{-%5.3f x} + %5.3f$' % tuple(popt), linewidth=2)
plt.xlabel('Time (ns)', fontsize=20)
plt.ylabel('SP', fontsize=20)
plt.title("Water Survival Probability, "+str(Model)+" $Ca^{+2}$ model, R"+str(Rep), fontsize=30)
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig("ST_"+str(Model)+"Ca_R"+str(Rep)+".jpg")

