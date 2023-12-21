from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


sci = fits.getdata('../data/lamp/20220409_lamp_stack_calibrated.fits')

# Get the 1-D spectrum
data = sci[1000:1200, :]   # Choose the central region
spectrum = np.mean(data, axis=0)

# Find peaks to identify the lines
peaks, _ = find_peaks(spectrum, height=10)
peaks = peaks[peaks <= 2700]   # Remove the last peak
pixel_coords = peaks.tolist()
print("Identified peak pixel coordinates:", pixel_coords)

known_wavelengths = [5852, 6562, 6965, 7067, 7272, 7383, 7503, 7514, 7635, 7723, 7948, 8006, 8103, 8115, 8264]
print(known_wavelengths)

for degree in range(11):
    coef = np.polyfit(pixel_coords, known_wavelengths, degree)
    poly_func = np.poly1d(coef)

    wavelengths = poly_func(pixel_coords)
    error = np.mean(np.abs(wavelengths - known_wavelengths))
    print(f"Degree {degree}: Mean absolute error = {error}")

