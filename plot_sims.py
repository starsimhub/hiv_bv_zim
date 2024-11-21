# %% Imports and settings
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
from utils import set_font, get_y

location = 'zimbabwe'


def plot_hiv_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='hiv_plots'):
    """ Create quantile or individual plots of HIV epi dynamics """
    set_font(size=20)
    fig, axes = pl.subplots(2, 2, figsize=(8, 7))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    hiv_data = pd.read_csv(f'data/{location}_hiv_data.csv')
    hiv_data = hiv_data.loc[(hiv_data.year >= start_year) & (hiv_data.year <= end_year)]

    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]

    # HIV infections
    pn = 0
    ax = axes[pn]
    resname = 'hiv.new_infections'
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='HIV infections')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV deaths
    ax = axes[pn]
    resname = 'hiv.new_deaths'
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='HIV deaths')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV-related deaths')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # PLHIV
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['hiv.n_infected'], color='k')  # label='UNAIDS',
    resnames = {'Total': 'hiv.n_infected', 'Dx': 'hiv.n_diagnosed', 'Treated': 'hiv.n_on_art'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = get_y(dfplot, which, rname)
        line, = ax.plot(x, y, label=rlabel)
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('PLHIV')
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV prevalence
    ax = axes[pn]
    resname = 'hiv.prevalence'
    ax.scatter(hiv_data.year, hiv_data[resname] * 100, label='Data', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y*100, label='Prevalence')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV prevalence (%)')
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year) + "_" + which + ".png", dpi=100)

    return fig


def plot_bv_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='sti_plots', fext=None):
    """ Create quantile or individual sim plots of BV """
    set_font(size=30)
    fig, axes = pl.subplots(3, 4, figsize=(25, 12))
    axes = axes.ravel()
    if which == 'multi': alphas = np.linspace(0.2, 0.5, len(percentile_pairs))
    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]
    pn = 0

    # Incidence
    ax = axes[pn]
    resnames = {'Total': 'bv.new_female_infections', 'Symptomatic': 'bv.new_female_symptomatic'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = get_y(dfplot, which, rname)
        line, = ax.plot(x, y, label=rlabel)
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())

    ax.set_title('BV incidence')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Burden
    ax = axes[pn]
    resnames = {'Total': 'bv.n_female_infected', 'Symptomatic': 'bv.n_female_symptomatic'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = get_y(dfplot, which, rname)
        line, = ax.plot(x, y, label=rlabel)
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())

    ax.set_title('BV burden')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Prevalence
    ax = axes[pn]
    resnames = {'Total': 'bv.female_adult_prevalence', 'Symptomatic': 'bv.female_symp_adult_prevalence'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = get_y(dfplot, which, rname)
        line, = ax.plot(x, y*100, label=rlabel)
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x, yl* 100, yu* 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('BV prevalence (%)')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year) + "_" + which + fext + ".png", dpi=100)

    return fig


if __name__ == '__main__':

    plot_single = False
    plot_multi = True

    if plot_multi:
        df_stats = sc.loadobj('results/multi_res_stats.df')
        percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]
        plot_hiv_sims(df_stats, start_year=2000, percentile_pairs=percentile_pairs)
        plot_bv_sims(df_stats, start_year=2000, percentile_pairs=percentile_pairs, which='multi')

    if plot_single:
        df = sc.loadobj('results/sim.df')
        plot_bv_sims(df, start_year=2000, which='single')
