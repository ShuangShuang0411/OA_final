from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


sci = fits.getdata('../data/M94/20220409_M94_stack_final.fits')

# Choose the area to plot the 1D spectrum
#data = sci[353:360, :]   # core 
data = sci[477:482, :]   # upper arm
#data = sci[215:220, :]   # lower arm
spectrum = np.mean(data, axis=0)

# The wavelength conversion function
def poly_func(x):
    return 9.533e-09*x**3 - 4.115e-05*x**2 + 1.88*x + 3347
wavelengths = poly_func(np.arange(len(spectrum)))

plt.plot(wavelengths, spectrum, linewidth=0.5)
plt.yscale('log')
plt.ylim(0.001, 6)
plt.xlabel('Wavelength (angstrom)')
plt.ylabel('Intensity')
plt.title('M94 upper arm spectrum (stack)')
plt.savefig("../figs_M94/spectrum_M94_stack_uarm.png", dpi=200)
plt.show()

