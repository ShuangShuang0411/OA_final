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
dark_files = sorted(glob('../data/dark/20220409_Dark*_dark600s.fits'))
master_dark = get_average_images(dark_files)
master_dark -= master_bias
fits.writeto('../data/dark/master_dark600s.fits', master_dark, overwrite=True)

# Flat
flat_files = sorted(glob('../data/flat/20220409_flat*.fits'))
master_flat = get_average_images(flat_files)
master_flat -= master_bias
master_flat = master_flat / np.nanmean(master_flat)
fits.writeto('../data/flat/master_flat.fits', master_flat, overwrite=True)

# Calibrate the science data
sci_files = sorted(glob('../data/M94position2/20220409_M94position2-*bin1_600s.fits'))
sci = get_average_images(sci_files)
hdu = fits.open('../data/M94position2/20220409_M94position2-001bin1_600s.fits')
hdr = hdu[0].header
exptime = hdr['EXPTIME']
#plot_image(sci)
#plt.savefig("../figs_M94position2/original_stack.png", dpi=200)

calibrated_sci = ((sci - master_bias - master_dark) / exptime) / master_flat 
#plot_image(calibrated_sci)
#plt.savefig("../figs_M94position2/calibrated_stack.png", dpi=200)

cut_sci = calibrated_sci[775:1350,:]
fits.writeto('../data/M94position2/20220409_M94position2_stack_calibrated.fits', cut_sci, overwrite=True)


# Calculate SNR
std = np.nanstd(cut_sci[400:540, 1500:2000])
print(std)

# Background subtraction
start_row = 400
end_row = 540
rect = plt.Rectangle((0, start_row), cut_sci.shape[1], end_row-start_row,
                      linewidth=1, edgecolor='red', facecolor='none')
plot_image(cut_sci)
plt.gca().add_patch(rect)
plt.savefig("../figs_M94position2/cut_bg_stack.png", dpi=200)
plt.show()
bg_reg = cut_sci[start_row:end_row, :]
bg = np.mean(bg_reg, axis=0)
cut_sci -= bg.reshape(1, -1)

# Rotate the image
point1 = (1211, 185)
point2 = (2325, 167)
plot_image(cut_sci)
plt.plot([point1[0], point2[0]], [point1[1], point2[1]], color='red')
plt.savefig("../figs_M94position2/cut_nobg_stack.png", dpi=200)
plt.show()
delta_y = point2[1] - point1[1]
delta_x = point2[0] - point1[0]
angle_radians = np.arctan2(delta_y, delta_x)
angle_degrees = np.degrees(angle_radians)

rotated_sci = ndimage.rotate(cut_sci, angle_degrees)
plot_image(rotated_sci)
plt.savefig("../figs_M94position2/cut_rot_stack.png", dpi=200)
plt.show()

fits.writeto('../data/M94position2/20220409_M94position2_stack_final.fits', rotated_sci, overwrite=True)

