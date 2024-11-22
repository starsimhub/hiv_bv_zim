"""
Create HIV model and interventions
"""

# %% Imports and settings
import numpy as np
import starsim as ss
import pandas as pd
import stisim as sti


def get_testing_products():
    """
    Define HIV products and testing interventions
    """
    years = np.arange(1990, 2021)  # Years for testing
    n_years = len(years)
    fsw_prob = np.linspace(0, 0.75, n_years)
    low_cd4_prob = np.linspace(0, 0.85, n_years)
    gp_prob = np.linspace(0, 0.5, n_years)

    # FSW agents who haven't been diagnosed or treated yet
    def fsw_eligibility(sim):
        return sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    fsw_testing = sti.HIVTest(
        years=years,
        test_prob_data=fsw_prob,
        name='fsw_testing',
        eligibility=fsw_eligibility,
        label='fsw_testing',
    )

    # Non-FSW agents who haven't been diagnosed or treated yet
    def other_eligibility(sim):
        return ~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    other_testing = sti.HIVTest(
        years=years,
        test_prob_data=gp_prob,
        name='other_testing',
        eligibility=other_eligibility,
        label='other_testing',
    )

    # Agents whose CD4 count is below 200.
    def low_cd4_eligibility(sim):
        return (sim.diseases.hiv.cd4 < 200) & ~sim.diseases.hiv.diagnosed

    low_cd4_testing = sti.HIVTest(
        years=years,
        test_prob_data=low_cd4_prob,
        name='low_cd4_testing',
        eligibility=low_cd4_eligibility,
        label='low_cd4_testing',
    )

    return fsw_testing, other_testing, low_cd4_testing


def make_hiv():
    """ Make HIV arguments for sim"""
    hiv = sti.HIV(
        beta={'structuredsexual': [1, 1], 'maternal': [1, 0.]},
        beta_m2f=0.006,
        beta_m2c=0.01,
        eff_condom=0.95,
        dur_on_art=ss.lognorm_ex(25, 5),
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=0.6,
    )
    return hiv


def make_hiv_intvs():

    n_art = pd.read_csv(f'data/n_art.csv').set_index('year')
    n_vmmc = pd.read_csv(f'data/n_vmmc.csv').set_index('year')
    fsw_testing, other_testing, low_cd4_testing = get_testing_products()
    art = sti.ART(coverage_data=n_art)
    vmmc = sti.VMMC(coverage_data=n_vmmc)

    interventions = [
        fsw_testing,
        other_testing,
        low_cd4_testing,
        art,
        vmmc
    ]

    return interventions
