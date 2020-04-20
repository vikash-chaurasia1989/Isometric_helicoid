import numpy as np
import math
import matplotlib
import scipy.optimize as optimize
# shortcuts
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import pi as PI
from scipy.optimize import fsolve
from scipy.optimize import minimize
from mayavi.mlab import *
#import plotly.graph_objects as go
from plotly.offline import init_notebook_mode, iplot
#import plotly.plotly as py
import plotly.graph_objs as go





def plotplytest(x,y,z):
    fig = go.Figure(data=[go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode= 'markers',
    marker=dict(
        size=2,
        color=z,                # set color to an array/list of desired values
        colorscale= 'Viridis',   # choose a colorscale
        opacity=1

        ),

        #colorbar=dict(thickness=10, tickvals=[np.min(z), np.max(z)], ticktext=['Low', 'High']),
        )])

    # tight layout
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    #fig.add_trace(colorbar_trace)
    fig.show()
    return
