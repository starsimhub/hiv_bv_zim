"""
Create a model with HIV plus 4 co-circulating discharging STIs:
    - chlamydia, gonorrhea, trichomoniasis, and other (BV+)
Used for evaluation of etiological tests compared to syndromic management.
"""

# %% Imports and settings
import sciris as sc
import starsim as ss
import stisim as sti
import pandas as pd

from hiv_model import make_hiv, make_hiv_intvs
from utils import unneeded_results


def make_bv(beta_m2f=0.2):
    bv = sti.BV(
        beta_m2f=beta_m2f,
        init_prev_data=pd.read_csv('data/init_prev_bv.csv'),
    )
    return bv


def make_sim(seed=1, n_agents=None, beta_m2f=0.15, dt=1/12, start=1980, stop=2030, debug=False, verbose=1/12):

    total_pop = {1970: 5.203e6, 1980: 7.05e6, 1990: 9980999, 2000: 11.83e6}[start]
    if n_agents is None: n_agents = [int(5e3), int(5e2)][debug]
    if dt is None: dt = [1/12, 1][debug]

    ####################################################################################################################
    # Demographic modules
    ####################################################################################################################
    fertility_rates = {'fertility_rate': pd.read_csv(f'data/asfr.csv')}
    pregnancy = ss.Pregnancy(pars=fertility_rates)
    death_rates = {'death_rate': pd.read_csv(f'data/deaths.csv')}
    death = ss.Deaths(death_rates)

    ####################################################################################################################
    # People and networks
    ####################################################################################################################
    ppl = ss.People(n_agents, age_data=pd.read_csv(f'data/age_dist_{start}.csv', index_col='age')['value'])
    sexual = sti.FastStructuredSexual(
        acts=ss.lognorm_ex(80, 30),
        prop_f1=0.2,
        prop_f2=0.05,
        prop_m1=0.2,
        f1_conc=0.05,
        f2_conc=0.25,
        m1_conc=0.15,
        m2_conc=0.3,
        p_pair_form=0.8,  # 0.6,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
    )
    maternal = ss.MaternalNet()

    ####################################################################################################################
    # Diseases
    ####################################################################################################################
    bv = make_bv(beta_m2f=beta_m2f)
    hiv = make_hiv()
    diseases = [bv, hiv]

    ####################################################################################################################
    # Interventions and analyzers
    ####################################################################################################################
    intvs = make_hiv_intvs()

    sim = ss.Sim(
        dt=dt,
        rand_seed=seed,
        total_pop=total_pop,
        start=start,
        stop=stop,
        people=ppl,
        diseases=diseases,
        networks=[sexual, maternal],
        demographics=[pregnancy, death],
        interventions=intvs,
        analyzers=[],
        connectors=sti.hiv_bv(hiv, bv),
        verbose=verbose,
    )

    return sim


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    do_run = True
    do_save = True

    if do_run:
        sim = make_sim(seed=seed, debug=debug, start=1980, stop=2030)
        sim.run(verbose=1/12)
        df = sti.finalize_results(sim, modules_to_drop=unneeded_results)
        if do_save: sc.saveobj('results/sim.df', df)

    # Process and plot
    from plot_sims import *
    df = sc.loadobj('results/sim.df')
    plot_sti_sims(df, start_year=2000, end_year=2040, which='single', fext='_alt')
    plot_hiv_sims(df, start_year=2000, which='single')

    print('Done.')


