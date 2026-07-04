import numpy as np
import pandas as pd
from scipy.optimize import least_squares

data = pd.read_csv("xy_data.csv")
xs, ys = data.iloc[:,0].values, data.iloc[:,1].values

def residuals(params):
    theta_deg, M, X = params
    theta = np.radians(theta_deg)

    u = (xs - X) * np.cos(theta) + (ys - 42) * np.sin(theta)
    v = -(xs - X) * np.sin(theta) + (ys - 42) * np.cos(theta)

    predicted_v = np.exp(M * np.abs(u)) * np.sin(0.3 * u)
    return v - predicted_v

res = least_squares(
    residuals,
    x0=[25, 0.0, 50],
    bounds=([0, -0.05, 0], [50, 0.05, 100]),
    xtol=1e-15, ftol=1e-15, gtol=1e-15  
)

theta_deg, M, X = res.x
print("theta (deg):", theta_deg)
print("theta (rad):", np.radians(theta_deg))
print("M:", M)
print("X:", X)
print("residual norm:", np.sum(res.fun**2))