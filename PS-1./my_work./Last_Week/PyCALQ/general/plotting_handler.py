#!/usr/bin/env python3
"""
Minimal plotting handler - without pylatex dependency
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
import numpy as np
import logging
from matplotlib.lines import Line2D

# Basic matplotlib settings
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.shadow'] = False
plt.rcParams['legend.framealpha'] = '0.8'
plt.rcParams['legend.facecolor'] = 'white'
plt.rcParams['legend.edgecolor'] = '0.5'
plt.rcParams['legend.fontsize'] = 15.0

# Try to use latex if available
try:
    plt.rcParams['text.usetex'] = True
except:
    plt.rcParams['text.usetex'] = False

class PlottingHandler:
    """Minimal plotting handler"""
    
    latex = plt.rcParams['text.usetex']
    
    def __init__(self):
        self.fig = None
        self.figwidth = 10
        self.figheight = 6.132
    
    def create_fig(self, figwidth, figheight):
        """Create a new figure"""
        self.fig = plt.figure(figsize=(figwidth, figheight))
        self.figheight = figheight
        self.figwidth = figwidth
        return self.fig
    
    def clf(self):
        """Clear the current figure"""
        plt.clf()
        if self.fig:
            self.fig.set_size_inches(self.figwidth, self.figheight)
    
    def set_figsize(self, figwidth, figheight):
        """Set figure size"""
        self.figheight = figheight
        self.figwidth = figwidth
        if self.fig:
            self.fig.set_size_inches(figwidth, figheight)
    
    def save_pdf(self, filename, transparent=True):
        """Save as PDF"""
        plt.savefig(filename, transparent=transparent, bbox_inches='tight')
    
    def save_pickle(self, filename):
        """Save as pickle (placeholder)"""
        logging.warning("Pickle saving not implemented in minimal version")
    
    def legend(self):
        """Show legend"""
        plt.legend()
    
    def ylim(self, ymin=None, ymax=None):
        """Get or set y limits"""
        if ymin is None and ymax is None:
            return plt.ylim()
        else:
            return plt.ylim((ymin, ymax))
    
    def set_y_logscale(self):
        plt.yscale('log')
    
    def set_y_linscale(self):
        plt.yscale('linear')
    
    def single_channel_plot(self, fig_params, channel, irreps, x_in, y_in):
        """Single channel scattering plot"""
        figwidth, figheight = fig_params
        plt.figure(figsize=(figwidth, figheight))
        
        markers = ['o', 's', '^', 'D', 'v', '*']
        legend_handles = []
        x, x_range = x_in
        y, y_range = y_in
        
        shape_count = 0
        for psq in irreps:
            for irrep in irreps[psq][0]:
                marker = markers[shape_count % len(markers)]
                shape_count += 1
                color = plt.cm.tab10(shape_count % 10)
                
                for level in irreps[psq][0][irrep]:
                    label = f"{irrep}({psq})" if len(irrep) > 1 else f"{irrep}"
                    plt.plot(x_range[psq][irrep][level], y_range[psq][irrep][level], 
                            color=color, alpha=0.5)
                    plt.plot(x[psq][irrep][level], y[psq][irrep][level], 
                            marker=marker, color=color, label=label)
                    legend_handles.append(Line2D([0], [0], marker=marker, 
                        color=color, markerfacecolor=color, markersize=10, label=label))
        
        plt.axhline(y=0, color='black')
        plt.axvline(x=0, color='black')
        plt.xlabel("$q^{*2} / m_{\pi}^2$", fontsize=14)
        plt.ylabel("$q^{*} / m_{\pi} \\cot \\delta $", fontsize=14)
        plt.legend(handles=legend_handles, loc='upper left', fontsize=10)
        plt.tight_layout()
    
    def summary_plot(self, indexes, levels, errs, xticks, reference=None, 
                     thresholds=[], label=None, index=0, ndatasets=1, shift=False, filled=True):
        """Spectrum summary plot"""
        plt.figure()
        plt.errorbar(x=indexes, y=levels, yerr=errs, fmt='o', capsize=5, label=label)
        plt.xticks(range(len(xticks)), [str(x) for x in xticks])
        plt.ylabel("Energy")
        plt.legend()
        plt.tight_layout()

# Make spectrum plotting function available
def make_spectrum_plot(*args, **kwargs):
    """Wrapper for spectrum plot"""
    handler = PlottingHandler()
    if 'energy_cm_data' in kwargs:
        # Simple spectrum plot
        data = kwargs.get('energy_cm_data', {})
        plt.figure(figsize=(10, 6))
        for psq, irreps in data.items():
            for irrep, levels in irreps.items():
                for level, energies in levels.items():
                    if isinstance(energies, list) and len(energies) > 0:
                        e_mean = np.mean(energies)
                        e_std = np.std(energies) if len(energies) > 1 else 0.01
                        plt.errorbar(f"{irrep}({psq})", e_mean, yerr=e_std, fmt='o', capsize=5)
        plt.ylabel("E_cm / m_N")
        plt.xticks(rotation=45)
        plt.tight_layout()
    return plt.gcf()

# Aliases for compatibility
apply_plot_style = lambda use_tex=True: None

__all__ = [
    'PlottingHandler',
    'make_spectrum_plot',
    'apply_plot_style',
]
