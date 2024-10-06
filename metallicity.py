from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Define file paths
data_file = "/Users/justinhew/Downloads/nihao_uhd_simulation_g8.26e11_xyz_positions_and_oxygen_ao.fits"
figures_dir = "./figures"

# Create directories if not exist
os.makedirs(figures_dir, exist_ok=True)
# with fits.open(data_file) as hdul:
#     data = hdul[1].data
#     print(data.columns)

# Load data from FITS file
with fits.open(data_file) as hdul:
    data = hdul[1].data
    x, y, z = data['x'], data['y'], data['z']
    ao = data['A_O']

# Calculate galactocentric radius
RGal = np.sqrt(x**2 + y**2 + z**2)

# Define a linear function for fitting
def linear_fit(x, a, b):
    return a * x + b

# Fit the data
popt, pcov = curve_fit(linear_fit, RGal, ao)
slope, intercept = popt
slope_uncertainty, intercept_uncertainty = np.sqrt(np.diag(pcov))

# Calculate the linear fit and residuals
fitted_ao = linear_fit(RGal, *popt)
residuals = ao - fitted_ao

# Plot 2-panel figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# (a) Logarithmic density plot of RGal vs A(O), with linear fit and legend
hb = ax1.hexbin(RGal, ao, gridsize=50, cmap='inferno', bins='log')
ax1.plot(RGal, fitted_ao, color='cyan', linewidth=2, label=f'Fit: A(O) = {slope:.3f} * RGal + {intercept:.3f}')
ax1.set_xlabel('Galactocentric Radius RGal (kpc)')
ax1.set_ylabel('Gas Phase Metallicity A(O)')
ax1.set_title('Logarithmic Density Plot of RGal vs A(O)')
ax1.legend()
cb = fig.colorbar(hb, ax=ax1)
cb.set_label('log(N)')

# (b) Residuals of the fit, RGal vs ∆A(O)
ax2.scatter(RGal, residuals, s=1, alpha=0.5, color='blue')
ax2.axhline(0, color='red', linestyle='--')
ax2.set_xlabel('Galactocentric Radius RGal (kpc)')
ax2.set_ylabel('Residuals ∆A(O)')
ax2.set_title('Residuals of Linear Fit (RGal vs ∆A(O))')

# Save 2-panel figure
fig.savefig(os.path.join(figures_dir, 'radial_metallicity_fit_and_residuals.png'))
plt.close(fig)

# Statistical metrics: Root Mean Square Error (RMSE)
rmse = np.sqrt(np.mean(residuals**2))

# Generate 2D histograms for the x vs y plane
bins = 50  # Define number of bins for 2D histogram

# (a) 2D-histogram of the median simulated A(O)
hist_a, xedges, yedges = np.histogram2d(x, y, bins=bins, weights=ao, density=False)
count_a, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
median_ao = np.divide(hist_a, count_a, out=np.zeros_like(hist_a), where=count_a != 0)

# (b) 2D-histogram of the median fitted A(O)
hist_f, _, _ = np.histogram2d(x, y, bins=[xedges, yedges], weights=fitted_ao, density=False)
median_fitted_ao = np.divide(hist_f, count_a, out=np.zeros_like(hist_f), where=count_a != 0)

# (c) 2D-histogram of the median residuals ∆A(O)
hist_r, _, _ = np.histogram2d(x, y, bins=[xedges, yedges], weights=residuals, density=False)
median_residuals = np.divide(hist_r, count_a, out=np.zeros_like(hist_r), where=count_a != 0)

# Plot 3-panel figure for 2D histograms
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# (a) Median simulated A(O)
im1 = axes[0].imshow(median_ao.T, origin='lower', aspect='auto', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap='viridis')
axes[0].set_title('Median Simulated A(O)')
axes[0].set_xlabel('x (kpc)')
axes[0].set_ylabel('y (kpc)')
fig.colorbar(im1, ax=axes[0])

# (b) Median fitted A(O)
im2 = axes[1].imshow(median_fitted_ao.T, origin='lower', aspect='auto', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap='viridis')
axes[1].set_title('Median Fitted A(O)')
axes[1].set_xlabel('x (kpc)')
fig.colorbar(im2, ax=axes[1])

# (c) Median residuals ∆A(O)
im3 = axes[2].imshow(median_residuals.T, origin='lower', aspect='auto', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap='coolwarm')
axes[2].set_title('Median Residuals ∆A(O)')
axes[2].set_xlabel('x (kpc)')
fig.colorbar(im3, ax=axes[2])

# Save 3-panel figure
fig.savefig(os.path.join(figures_dir, '2D_histograms_of_metallicity.png'))
plt.close(fig)

# Output summary of fitting results
{
    "slope": slope,
    "slope_uncertainty": slope_uncertainty,
    "intercept": intercept,
    "intercept_uncertainty": intercept_uncertainty,
    "rmse": rmse
}


