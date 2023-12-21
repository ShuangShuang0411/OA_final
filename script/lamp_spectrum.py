from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


sci = fits.getdata('../data/lamp/20220409_lamp_stack_calibrated.fits')

# Get the 1-D spectrum
data = sci[1000:1200, :]   # Choose the central region
spectrum = np.mean(data, axis=0)

plt.plot(spectrum, linewidth=1.0)
plt.yscale('log')
plt.ylim(0.1, 600)
plt.xlabel('x (pixel)')
plt.ylabel('Intensity')
plt.title('Ne Ar lamp spectrum')
plt.savefig("../figs_lamp/spectrum_org.png", dpi=200)

# Find peaks to identify the lines
peaks, _ = find_peaks(spectrum, height=10)
peaks = peaks[peaks <= 2700]   # Remove the last peak
pixel_coords = peaks.tolist()
print("Identified peak pixel coordinates:", pixel_coords)

# Check the chosen lines
#peak_spectrum = np.zeros_like(spectrum)
#peak_spectrum[pixel_coords] = 1
#plt.clf()
#plt.plot(peak_spectrum, linewidth=1.0)
#plt.savefig("chosen_lines.png", dpi=200)

known_wavelengths = [5852, 6562, 6965, 7067, 7272, 7383, 7503, 7514, 7635, 7723, 7948, 8006, 8103, 8115, 8264]
print(known_wavelengths)

# Fit a polynomial to connent pixel coords to wavelengths
degree = 3 
coef = np.polyfit(pixel_coords, known_wavelengths, degree)
poly_func = np.poly1d(coef)
print("conversion fnt:", poly_func)

# Convert pixel coords to wavelengths using the polynomial function
wavelengths = poly_func(np.arange(len(spectrum)))

plt.clf()
plt.plot(wavelengths, spectrum, linewidth=1.0)
plt.yscale('log')
plt.ylim(0.1, 600)
plt.xlabel('Wavelength (angstrom)')
plt.ylabel('Intensity')
plt.title('Ne Ar lamp spectrum')
plt.savefig("../figs_lamp/spectrum_final.png", dpi=200)

