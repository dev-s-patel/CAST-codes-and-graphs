import matplotlib.pyplot as plt
from getdist import plots, loadMCSamples

# Load chains from the chains/ directory
samples = loadMCSamples('./chains/cast_planck', settings={'ignore_rows': 0.3})

# Publication styling
g = plots.get_subplot_plotter(width_inch=8.5)
g.settings.figure_legend_frame = False
g.settings.alpha_filled_add = 0.6
g.settings.axes_fontsize = 11
g.settings.lab_fontsize = 12

# Parameters to display
params_to_plot = ['xi_cast', 'k_c', 'H0', 'ombh2', 'omch2']

# Triangle plot with 1-sigma (68%) & 2-sigma (95%) contours + theoretical marker
g.triangle_plot(
    samples, 
    params_to_plot, 
    filled=True, 
    contour_colors=['#1f77b4'],
    diag1d_kwargs={'color': '#1f77b4', 'lw': 1.8},
    markers={'xi_cast': 0.04955}
)

# Save figures
plt.savefig('cast_planck_contours.pdf', bbox_inches='tight', dpi=300)
plt.savefig('cast_planck_contours.png', bbox_inches='tight', dpi=300)
print("Plots saved successfully as cast_planck_contours.pdf and cast_planck_contours.png")
