'''
Simulation of beam interaction with coupling impedance and damper
for a single bunch
'''

# using Python 2.7:
from __future__ import division, print_function
range_ = range
range = xrange

import sys
# PyHEADTAIL location if it's not already in the PYTHONPATH environment variable
# sys.path.append('/home/XYZ')



import time

import numpy as np
np.random.seed(10000042)
import h5py
from scipy.constants import e, m_p, c

from PyHEADTAIL.particles.slicing import UniformBinSlicer
from PyHEADTAIL.impedances.wakes import WakeTable, WakeField
from PyHEADTAIL.feedback.transverse_damper import TransverseDamper
from PyHEADTAIL.monitors.monitors import (
    BunchMonitor, ParticleMonitor, SliceMonitor)


n_macroparticles = 1000000 # number of macro-particles to resolve the beam
n_turns = 600000 # simulation time
n_turns_slicemon = 2048 # recording span of the slice statistics monitor


# COMMENT THE NON-WANTED SET-UP:

# # injection
# machine_configuration = 'LHC-injection'
# wakefile = ('./wakeforhdtl_PyZbase_Allthemachine_450GeV'
#             '_B1_LHC_inj_450GeV_B1.dat')

# flat-top
machine_configuration = 'LHC_6.5TeV_collision_2016'
wakefile = ('./wakeforhdtl_PyZbase_Allthemachine_6p5TeV'
            '_B1_LHC_ft_6.5TeV_B1.dat')

# ---

def get_nonlinear_params(chroma, i_oct, p0=6.5e12*e/c):
    '''Arguments:
        - chroma: first-order chromaticity Q'_{x,y}, identical
          for both transverse planes
        - i_oct: octupole current in A (positive i_oct means
          LOF = i_oct > 0 and LOD = -i_oct < 0)
    '''
    # factor 2p0 is PyHEADTAIL's convention for d/dJx instead of
    # MAD-X's convention of d/d(2Jx)
    app_x = 2 * p0 * 27380.10941 * i_oct / 100.
    app_y = 2 * p0 * 28875.03442 * i_oct / 100.
    app_xy = 2 * p0 * -21766.48714 * i_oct / 100.
    Qpp_x = 4889.00298 * i_oct / 100.
    Qpp_y = -2323.147896 * i_oct / 100.
    return {
        'app_x': app_x,
        'app_y': app_y,
        'app_xy': app_xy,
        'Qp_x': [chroma,],# Qpp_x],
        'Qp_y': [chroma,],# Qpp_y],
        # second-order chroma commented out above!
    }


def run(intensity, chroma=0, i_oct=0):
    '''Arguments:
        - intensity: integer number of charges in beam
        - chroma: first-order chromaticity Q'_{x,y}, identical
          for both transverse planes
        - i_oct: octupole current in A (positive i_oct means
          LOF = i_oct > 0 and LOD = -i_oct < 0)
    '''


    # BEAM AND MACHINE PARAMETERS
    # ============================
    from LHC import LHC
    # energy set above will enter get_nonlinear_params p0
    assert machine_configuration == 'LHC_6.5TeV_collision_2016'
    machine = LHC(n_segments=1,
                  machine_configuration=machine_configuration,
                  **get_nonlinear_params(chroma=chroma, i_oct=i_oct))


    # BEAM
    # ====
    epsn_x = 3.e-6 # normalised horizontal emittance
    epsn_y = 3.e-6 # normalised vertical emittance
    sigma_z = 1.2e-9 * machine.beta*c/4. # RMS bunch length in meters

    bunch = machine.generate_6D_Gaussian_bunch_matched(
        n_macroparticles, intensity, epsn_x, epsn_y, sigma_z=sigma_z)

    print ("\n--> Bunch length and emittance: {:g} m, {:g} eVs.".format(
            bunch.sigma_z(), bunch.epsn_z()))


    # CREATE BEAM SLICERS
    # ===================
    slicer_for_slicemonitor = UniformBinSlicer(
        50, z_cuts=(-3*sigma_z, 3*sigma_z))
    slicer_for_wakefields = UniformBinSlicer(
        500, z_cuts=(-3*sigma_z, 3*sigma_z))


    # CREATE WAKES
    # ============
    wake_table1 = WakeTable(wakefile,
                            ['time', 'dipole_x', 'dipole_y',
                             'quadrupole_x', 'quadrupole_y',
                             # 'noquadrupole_x', 'noquadrupole_y',
                             'dipole_xy', 'dipole_yx',
                             # 'nodipole_xy', 'nodipole_yx',
                            ])
    wake_field = WakeField(slicer_for_wakefields, wake_table1)


    # CREATE DAMPER
    # =============
    dampingrate = 50
    damper = TransverseDamper(dampingrate, dampingrate)


    # CREATE MONITORS
    # ===============
    try:
        bucket = machine.longitudinal_map.get_bucket(bunch)
    except AttributeError:
        bucket = machine.rfbucket

    simulation_parameters_dict = {
        'gamma'    : machine.gamma,
        'intensity': intensity,
        'Qx'       : machine.Q_x,
        'Qy'       : machine.Q_y,
        'Qs'       : bucket.Q_s,
        'beta_x'   : bunch.beta_Twiss_x(),
        'beta_y'   : bunch.beta_Twiss_y(),
        'beta_z'   : bucket.beta_z,
        'epsn_x'   : bunch.epsn_x(),
        'epsn_y'   : bunch.epsn_y(),
        'sigma_z'  : bunch.sigma_z(),
    }
    bunchmonitor = BunchMonitor(
        outputpath+'/bunchmonitor_{:04d}_chroma={:g}'.format(it, chroma),
        n_turns, simulation_parameters_dict,
        write_buffer_to_file_every=512,
        buffer_size=4096)
    slicemonitor = SliceMonitor(
        outputpath+'/slicemonitor_{:04d}_chroma={:g}'.format(it, chroma),
        n_turns_slicemon,
        slicer_for_slicemonitor, simulation_parameters_dict,
        write_buffer_to_file_every=1, buffer_size=n_turns_slicemon)


    # TRACKING LOOP
    # =============
    machine.one_turn_map.append(damper)
    machine.one_turn_map.append(wake_field)


    # for slice statistics monitoring:
    s_cnt = 0
    monitorswitch = False

    print ('\n--> Begin tracking...\n')

    # GO!!!
    for i in range(n_turns):

        t0 = time.clock()

        # track the beam around the machine for one turn:
        machine.track(bunch)

        ex, ey, ez = bunch.epsn_x(), bunch.epsn_y(), bunch.epsn_z()
        mx, my, mz = bunch.mean_x(), bunch.mean_y(), bunch.mean_z()

        # monitor the bunch statistics (once per turn):
        bunchmonitor.dump(bunch)

        # if the centroid becomes unstable (>1cm motion)
        # then monitor the slice statistics:
        if not monitorswitch:
            if mx > 1e-2 or my > 1e-2 or i > n_turns - n_turns_slicemon:
                print ("--> Activate slice monitor")
                monitorswitch = True
        else:
            if s_cnt < n_turns_slicemon:
                slicemonitor.dump(bunch)
                s_cnt += 1

        # stop the tracking as soon as we have not-a-number values:
        if not all(np.isfinite(c) for c in [ex, ey, ez, mx, my, mz]):
            print ('*** STOPPING SIMULATION: non-finite bunch stats!')
            break

        # print status all 1000 turns:
        if i % 1000 == 0:
            t1 = time.clock()
            print ('Emittances: ({:.3g}, {:.3g}, {:.3g}) '
                   '& Centroids: ({:.3g}, {:.3g}, {:.3g})'
                   '@ turn {:d}, {:g} ms, {:s}'.format(
                        ex, ey, ez, mx, my, mz, i, (t1-t0)*1e3, time.strftime(
                            "%d/%m/%Y %H:%M:%S", time.localtime()))
            )

    print ('\n*** Successfully completed!')


if __name__ == '__main__':
    # iteration, attached to monitor name:
    it = 0
    # outputpath relative to this file:
    outputpath = './'

    # run the simulation:
    run(intensity=1.1e11, chroma=12, i_oct=0)
