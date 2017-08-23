"""Quick-and-dirty plots of :func:`~fluxpart.fluxpart.flux_partition` output"""

import matplotlib.pyplot as plt


PLTSTYLE = 'bmh'


def plot_hfdata(hfdata):
    """Note: converts q, c, T, and P from SI to common units"""
    plt.style.use(PLTSTYLE)
    fig, ax = plt.subplots(7, 1, figsize=(5, 8), sharex='all')
    ax[0] = _tsplot(hfdata['u'], r'u [m/s]', ax[0])
    ax[1] = _tsplot(hfdata['v'], r'v [m/s]', ax[1])
    ax[2] = _tsplot(hfdata['w'], r'w [m/s]', ax[2])
    ax[3] = _tsplot(1e3 * hfdata['q'], r'q [g/m^3]', ax[3])
    ax[4] = _tsplot(1e6 * hfdata['c'], r'c [mg/m^3]', ax[4])
    ax[5] = _tsplot(hfdata['T'] - 273.15, r'T [C]', ax[5])
    ax[6] = _tsplot(1e-3 * hfdata['P'], r'P [kPa]', ax[6])
    fig.tight_layout(pad=1, w_pad=0, h_pad=0)
    return fig


def _tsplot(series, yaxlab, ax=None):
    if ax is None:
        ax = plt.gca()
    ax.plot(series, linewidth=0.5)
    ax.set_ylabel('{} $\mathrm{{ {} }}$'.format(*yaxlab.split(maxsplit=1)))
    ax.set_xticklabels([])
    ax.locator_params(axis='y', tight=True, nbins=4)
    return ax


def plot_fluxes(fp_results, figtitle=None):
    """Note: converts Fq and Fc from SI to common units"""
    plt.style.use(PLTSTYLE)
    fcr = [1e6 * r['fluxes'].Fcr_mol for r in fp_results]
    fcp = [1e6 * r['fluxes'].Fcp_mol for r in fp_results]
    fc = [1e6 * r['fluxes'].Fc_mol for r in fp_results]
    fqe = [1e3 * r['fluxes'].Fqe_mol for r in fp_results]
    fqt = [1e3 * r['fluxes'].Fqt_mol for r in fp_results]
    fq = [1e3 * r['fluxes'].Fq_mol for r in fp_results]
    times = [r['label'] or i for i, r in enumerate(fp_results)]
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(6, 6), sharex='all')
    ax0.plot(times, fc, label='Total')
    ax0.plot(times, fcr, 'o', label='Respiration')
    ax0.plot(times, fcp, '^', label='Photosynthesis')
    ax1.plot(times, fq, label='Total')
    ax1.plot(times, fqe, 'o', label='Evaporation')
    ax1.plot(times, fqt, '^', label='Transpiration')
    ax0.set_ylabel(r'CO2 $[\mathrm{umol/m^2/s}]$')
    ax1.set_ylabel(r'H2O $[\mathrm{mmol/m^2/s}]$')
    ax0.legend(loc='best', fontsize=10)
    ax1.legend(loc='best', fontsize=10)
    if figtitle:
        fig.suptitle(figtitle)
    plt.gcf().autofmt_xdate()
    fig.tight_layout()
    return fig


def plot_fluxes_mass(fp_results, figtitle=None):
    """Note: converts Fq and Fc from SI to common units"""
    plt.style.use(PLTSTYLE)
    fcr = [1e6 * r['fluxes'].Fcr for r in fp_results]
    fcp = [1e6 * r['fluxes'].Fcp for r in fp_results]
    fc = [1e6 * r['fluxes'].Fc for r in fp_results]
    fqe = [1e3 * r['fluxes'].Fqe for r in fp_results]
    fqt = [1e3 * r['fluxes'].Fqt for r in fp_results]
    fq = [1e3 * r['fluxes'].Fq for r in fp_results]
    times = [r['label'] or i for i, r in enumerate(fp_results)]
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(6, 6), sharex='all')
    ax0.plot(times, fc, label='Total')
    ax0.plot(times, fcr, 'o', label='Respiration')
    ax0.plot(times, fcp, '^', label='Photosynthesis')
    ax1.plot(times, fq, label='Total')
    ax1.plot(times, fqe, 'o', label='Evaporation')
    ax1.plot(times, fqt, '^', label='Transpiration')
    ax0.set_ylabel(r'CO2 $[\mathrm{mg/m^2/s}]$')
    ax1.set_ylabel(r'H2O $[\mathrm{g/m^2/s}]$')
    ax0.legend(loc='best', fontsize=10)
    ax1.legend(loc='best', fontsize=10)
    if figtitle:
        fig.suptitle(figtitle)
    plt.gcf().autofmt_xdate()
    fig.tight_layout()
    return fig
