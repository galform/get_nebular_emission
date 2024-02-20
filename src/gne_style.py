from cycler import cycler

# see http://matplotlib.org/users/customizing.html for all options

style1 = {
    # Line styles
    'lines.linewidth': 2.5,
    'lines.antialiased': True,

    # Font
    'font.size': 16.0,
    'font.family': 'sans-serif',

    # Axes
    'axes.linewidth': 1.5,
    'axes.titlesize': 'small',
    'axes.labelsize': 'medium',
    'axes.prop_cycle': cycler('color', [
        '#1f77b4',  # blue
        '#ff7f0e',  # orange
        '#2ca02c',  # green
        '#d62728',  # red
        '#9467bd',  # purple
        '#8c564b',  # brown
        '#e377c2',  # magenta
        '#7f7f7f',  # gray
        '#bcbd22',  # yellow
        '#17becf',  # cyan
    ]),

    # Ticks
    'xtick.major.size': 6,
    'xtick.minor.size': 3,
    'xtick.major.width': 1.5,
    'xtick.minor.width': 1.5,
    'xtick.major.pad': 6,
    'xtick.minor.pad': 3,
    'xtick.labelsize': 'medium',
    'xtick.direction': 'in',
    'xtick.top': True,
    'xtick.bottom': True,
    'xtick.minor.visible': True,

    'ytick.major.size': 6,
    'ytick.minor.size': 3,
    'ytick.major.width': 1.5,
    'ytick.minor.width': 1.5,
    'ytick.major.pad': 6,
    'ytick.minor.pad': 3,
    'ytick.labelsize': 'medium',
    'ytick.direction': 'in',
    'ytick.right': True,
    'ytick.left': True,
    'ytick.minor.visible': True,

    # Legend
    'legend.fancybox': False,
    'legend.fontsize': 'small',
    'legend.scatterpoints': 1,
    'legend.numpoints': 1,
    'legend.loc': 'best',

    # Figure
    'figure.figsize': [6.5, 7],
    'figure.titlesize': 'large',

    # Images
    'image.cmap': 'magma',
    'image.origin': 'lower',

    # Saving
    'savefig.bbox': 'tight',
    'savefig.format': 'pdf',
}

# https://matplotlib.org/users/customizing.html
