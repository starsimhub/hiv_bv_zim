import starsim as ss
import sciris as sc
import pandas as pd
from model import make_sim
from calibration import Calibration

# Settings
debug = False
do_save = True

# Run settings for calibration (dependent on debug)
n_trials = [1000, 1][debug]  # How many trials to run for calibration
n_workers = [50, 1][debug]  # How many cores to use
storage = ["mysql://hpvsim_user@localhost/hivsim_db", None][debug]  # Storage for calibrations


def build_sim(sim, calib_pars):
    """ Logic for converting calibration parameters into simulation parameters """

    hiv = sim.diseases.hiv
    bv = sim.diseases.bv

    for k, pars in calib_pars.items():  # Loop over the calibration parameters
        if k == 'rand_seed':
            sim.pars.rand_seed = v
            continue

        v = pars['value']
        if k == 'hiv_beta_m2f':
            hiv.pars.beta_m2f = v
        elif k == 'p_base':
            bv.pars.p_base = v
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

    return sim


def run_calib(calib_pars=None):
    sc.heading('Beginning calibration')

    # Make the sim and data
    sim = make_sim()
    data = pd.read_csv('data/zimbabwe_calib_data.csv')

    # Make the calibration
    calib = Calibration(
        sim=sim,
        data=data,
        calib_pars=calib_pars,
        build_fn = build_sim,
        total_trials = n_trials,
        n_workers = n_workers,
        die = True,
        reseed=False,
        debug = debug,
    )

    # Perform the calibration
    sc.printcyan('\nPeforming calibration...')
    calib.calibrate()
    return sim, calib


#%% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    # Define the calibration parameters
    calib_pars = dict(
        hiv_beta_m2f=dict(low=0.001, high=0.1, guess=0.04),
        p_base=dict(low=0.01, high=0.1, guess=0.05),
    )

    sim, calib = run_calib(calib_pars=calib_pars)
    from utils import shrink_calib
    cal = shrink_calib(calib, n_results=100)
    sc.saveobj('results/calib.obj', cal)

    sc.toc(T)
    print('Done.')


