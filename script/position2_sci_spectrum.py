from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


sci = fits.getdata('../data/M94position2/20220409_M94position2_stack_final.fits')

# Choose the area to plot the 1D spectrum
data = sci[202:207, :]   # core 
#data = sci[324:328, :]   # upper arm
#data = sci[72:77, :]   # lower arm
spectrum = np.mean(data, axis=0)

# The wavelength conversion function
def poly_func(x):
    return 9.533e-09*x**3 - 4.115e-05*x**2 + 1.88*x + 3347
wavelengths = poly_func(np.arange(len(spectrum)))

plt.plot(wavelengths, spectrum, linewidth=0.5)
plt.yscale('log')
plt.ylim(0.05, 6)
plt.xlabel('Wavelength (angstrom)')
plt.ylabel('Intensity')
plt.title('M94 core spectrum (stack)')
plt.savefig("../figs_M94position2/spectrum_M94_stack_core.png", dpi=200)
plt.show()


