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


for i in range(1, 5):
#    file_name = f"../data/M94/20220409_M94-{i:03d}bin1_600s.fits"
#    file_name = f"../data/M94position2/20220409_M94position2-{i:03d}bin1_600s.fits"
    file_name = f"../data/lamp/20220409_lamp-{i:03d}bin1_150s.fits"
    hdu = fits.open(file_name)
    sci = hdu[0].data
    hdr = hdu[0].header
    exptime = hdr['EXPTIME']
    plot_image(sci)
    plt.savefig(f"../figs_lamp/original{i}.png", dpi=200)

