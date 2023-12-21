from astropy.io import fits
from astropy.visualization import ZScaleInterval
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage


def get_average_images(files):
    data = [fits.getdata(file) for file in files]
    image = np.nanmean(data, axis=0)
    return image

def plot_image(data):
    plt.clf()
    vmin, vmax = ZScaleInterval().get_limits(data)
    plt.imshow(data, origin='lower', vmin=vmin, vmax=vmax)
#    plt.imshow(data, origin='lower', vmin=-0.1, vmax=0.25)
    plt.colorbar()


# Bias
bias_files = sorted(glob('../data/bias/20220409_bias*.fits'))
master_bias = get_average_images(bias_files) 
fits.writeto('../data/bias/master_bias.fits', master_bias, overwrite=True)

# Dark
dark_files = sorted(glob('../data/dark/20220409_dark*_dark150s.fits'))
master_dark = get_average_images(dark_files)
master_dark -= master_bias
fits.writeto('../data/dark/master_dark150s.fits', master_dark, overwrite=True)

# Flat
flat_files = sorted(glob('../data/flat/20220409_flat*.fits'))
master_flat = get_average_images(flat_files)
master_flat -= master_bias
master_flat = master_flat / np.nanmean(master_flat)
fits.writeto('../data/flat/master_flat.fits', master_flat, overwrite=True)

# Calibrate the science data
sci_files = sorted(glob('../data/lamp/20220409_lamp-*bin1_150s.fits'))
sci = get_average_images(sci_files)
hdu = fits.open('../data/lamp/20220409_lamp-001bin1_150s.fits')
hdr = hdu[0].header
exptime = hdr['EXPTIME']
plot_image(sci)
plt.savefig("../figs_lamp/original_stack.png", dpi=200)

calibrated_sci = ((sci - master_bias - master_dark) / exptime)   #/ master_flat 
plot_image(calibrated_sci)
plt.savefig("../figs_lamp/calibrated_stack.png", dpi=200)

fits.writeto('../data/lamp/20220409_lamp_stack_calibrated.fits', calibrated_sci, overwrite=True)

