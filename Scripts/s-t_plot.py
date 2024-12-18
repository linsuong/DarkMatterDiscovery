import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# Experimental values and uncertainties
S_exp, T_exp = 0.06, 0.1  # Central values
sigma_S, sigma_T = 0.09, 0.07  # Uncertainties
rho = 0.91  # Correlation coefficient

# Function to plot the experimental ellipse
def plot_ellipse(ax, S_exp, T_exp, sigma_S, sigma_T, rho, n_std=1, **kwargs):
    cov = np.array([[sigma_S**2, rho * sigma_S * sigma_T],
                    [rho * sigma_S * sigma_T, sigma_T**2]])
    vals, vecs = np.linalg.eigh(cov)
    width, height = 2 * n_std * np.sqrt(vals)  # Ellipse axes
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))  # Ellipse rotation
    ellipse = Ellipse(xy=(S_exp, T_exp), width=width, height=height, angle=angle, **kwargs)
    ax.add_patch(ellipse)

# Theoretical predictions (example points)
S_theory = [0.0, 0.05, 0.1, 0.2]  # Example S values
T_theory = [0.1, 0.15, 0.1, 0.05]  # Corresponding T values

# Plotting
fig, ax = plt.subplots(figsize=(8, 6))

# Plot experimental ellipse
plot_ellipse(ax, S_exp, T_exp, sigma_S, sigma_T, rho, n_std=1, edgecolor='red', facecolor='none', linewidth=2, label="1σ EWPT")

# Plot theoretical points
ax.scatter(S_theory, T_theory, color='blue', s=50, label="Model Predictions")

# Labels and legend
ax.axhline(0, color='gray', linestyle='--')  # y = 0 line
ax.axvline(0, color='gray', linestyle='--')  # x = 0 line
ax.set_xlabel("S Parameter")
ax.set_ylabel("T Parameter")
ax.set_title("S-T Plane with EWPT Constraints")
ax.legend()
plt.grid()
plt.show()