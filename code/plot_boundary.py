"""Render Figure 1: the attainable region and the checkerboard comparison (requires matplotlib).

Run from the package root: python3 code/plot_boundary.py
Reads results/verification.csv; writes figures/boundary-figure.pdf and .png.
"""
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_boundary import curve  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
INK = '#362f3a'
MULBERRY = '#883958'
OCHRE = '#ab6e24'
MOSS = '#507160'
FILL = '#e8e3d2'


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
                         'axes.edgecolor':INK,'text.color':INK,
                         'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,
                         'pdf.fonttype':42,'ps.fonttype':42})
    fig, (region, diagnostics) = plt.subplots(1,2,figsize=(9.3,3.6),
        gridspec_kw={'width_ratios':[1.35,1]},layout='constrained')
    slopes = np.geomspace(.001,128,180)
    points = np.array([curve(k,4096) for k in slopes])
    x = np.r_[0,points[:,0],.5,1]
    lower = np.r_[0,points[:,1],-.5,-.5]
    xi_uniform = np.linspace(0,1,1001)
    region.fill_between(xi_uniform,np.interp(xi_uniform,x,lower),
                        np.sqrt(xi_uniform),color=FILL)
    region.plot(xi_uniform,np.sqrt(xi_uniform),color=MOSS,lw=1.7)
    region.plot(x,lower,color=MULBERRY,lw=2)
    region.plot([0,1],[0,-.5],color=OCHRE,lw=1.2,ls='--')
    region.annotate(r'$\xi=-2\psi$ (SD copulas)',(.5,-.2),xytext=(0,0),
                    textcoords='offset points',color=OCHRE,fontsize=8)
    region.scatter([0,.5,1],[0,-.5,-.5],s=15,color=MULBERRY,zorder=4)
    region.annotate(r'$\psi=\sqrt{\xi}$',(.55,.74),xytext=(-12,7),
                    textcoords='offset points',color=MOSS)
    region.annotate(r'$k=1$',curve(1,4096),xytext=(-17,-20),
                    textcoords='offset points',color=MULBERRY,
                    arrowprops={'arrowstyle':'-','color':MULBERRY,'lw':.7})
    region.text(.47,.11,'Attainable pairs',ha='center')
    region.set(xlabel=r"Chatterjee’s $\xi$",ylabel=r"Spearman’s $\psi$",
               xlim=(-.015,1.015),ylim=(-.64,1.05),xticks=np.arange(0,1.01,.2),
               yticks=[-.5,-.25,0,.25,.5,.75,1])
    region.set_title('A   Attainable region',loc='left',fontsize=10,weight='bold')
    data = np.genfromtxt(ROOT/'results'/'verification.csv',names=True,delimiter=',')
    for k,color,mark in [(.125,MOSS,'o'),(1,MULBERRY,'s'),(8,OCHRE,'^')]:
        d = data[data['k']==k]
        diagnostics.loglog(d['n'],d['grid_minus_continuum_J'],marker=mark,
                           color=color,lw=1.6,ms=4,label=fr'$k={k:g}$')
    diagnostics.set(xlabel='Checkerboard cells per axis',
                    ylabel=r'Excess $\xi+2k\psi$ above boundary',
                    xticks=[16,32,64],xlim=(14,75))
    diagnostics.xaxis.set_major_formatter(ScalarFormatter())
    diagnostics.xaxis.set_minor_locator(plt.NullLocator())
    diagnostics.legend(frameon=False,fontsize=8,loc='upper right')
    diagnostics.set_title('B   Independent minimisation',loc='left',fontsize=10,weight='bold')
    for ax in [region,diagnostics]:
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y',color='#b7aea2',alpha=.28,lw=.6)
        ax.set_axisbelow(True)
    (ROOT/'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT/'figures'/'boundary-figure.pdf',bbox_inches='tight',metadata={'CreationDate':None})
    fig.savefig(ROOT/'figures'/'boundary-figure.png',dpi=220,bbox_inches='tight')


if __name__=='__main__':
    main()
