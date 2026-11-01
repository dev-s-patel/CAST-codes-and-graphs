import getdist.plots as gdplots
from getdist import loadMCSamples
import matplotlib.pyplot as plt

# Load Cobaya chains
samples = loadMCSamples('chains/cast_planck', settings={'ignore_rows': 0.1})

# Initialize triangle plot
g = gdplots.get_subplot_plotter(width_inch=8.5)
g.settings.axes_fontsize = 9
g.settings.lab_fontsize = 11
g.settings.legend_fontsize = 10

# Primary parameter subset
plot_params = ['ombh2', 'omch2', 'H0', 'ns', 'k_c', 'xi_cast']

g.triangle_plot(
    samples, 
    plot_params, 
    filled=True, 
    contour_colors=['#1b7837'],
    line_args={'lw': 1.6, 'color': '#1b7837'}
)

# Add theoretical prediction line on xi_cast 1D posterior
ax_xi = g.subplots[-1, -1]
ax_xi.axvline(0.050, color='crimson', linestyle='--', linewidth=1.5, label=r'$\xi_{\rm theory} \approx 0.050$')
ax_xi.legend(loc='upper right', frameon=False, fontsize=8)

plt.savefig('cast_planck_contours.pdf', bbox_inches='tight', dpi=300)
plt.savefig('cast_planck_contours.png', bbox_inches='tight', dpi=300)
print("[OK] Exported cast_planck_contours.pdf and .png")