# PyHEADTAIL-playground
A collection of jupyter notebooks and scripts for [PyHEADTAIL](https://github.com/PyCOMPLETE/PyHEADTAIL).

<a href="https://cern.ch/swanserver/cgi-bin/go?projurl=https://github.com/PyCOMPLETE/PyHEADTAIL-playground.git" target="_blank"><img src="http://swanserver.web.cern.ch/swanserver/images/badge_swan_white_150.png"/></a>

A quick start tutorial explaining the basics of `PyHEADTAIL` can be found on the top level of this repository:
[Quick Start Tutorial for PyHEADTAIL](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/Tutorial.ipynb)

## Full Simulation Examples

In `full_simulations/` you find example scripts for full `PyHEADTAIL` simulations, such as
* [LHC head-tail instability](https://github.com/PyCOMPLETE/PyHEADTAIL-playground/tree/master/full_simulations/LHC_instability)

## Simulation Notebooks

In `simulation_notebooks/` you find jupyter notebooks running and explaining concrete `PyHEADTAIL` simulation examples, such as
* [2-particle model explaining the transverse mode coupling instability](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/simulation_notebooks/Transverse_Mode_Coupling_Instability/TMCI_2particle_model.ipynb) with instructive animations
* [transverse ideal resistive and reactive damper](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/simulation_notebooks/TransverseDamper.ipynb) with some exponential rise time fitting
* [triple splitting as in the CERN Proton Synchrotron RF gymnastics](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/simulation_notebooks/PS-TripleBunchSplitting.ipynb)
* [space charge tutorial](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/simulation_notebooks/SpaceChargeTutorial.ipynb) with incoherent footprints in tune diagram

## How To Notebooks

With a bit more detail and focus, modules and concepts of `PyHEADTAIL` are explained in jupyter notebooks in `howto_notebooks/`. These how-to's deal with
* [aperture and loss modelling](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/ApertureNLossesTest.ipynb)
* [optics function values inferred from beam distribution](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/BeamOptics.ipynb)
* [detuning models for chromaticity and octupole amplitude detuning](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/DetunersTest.ipynb)
* [generation of macro-particles in different ways](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/GeneratorTest.ipynb)
* [monitoring of beam distribution during simulation](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/MonitorTest.ipynb)
* [thin multipole kicks](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/MultipolesTest.ipynb)
* [PyHEADTAIL on the GPU tutorial](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/PyHEADTAIL_on_GPU_Tutorial.ipynb)
* [matching of arbitrary stationary distributions into arbitrary RF buckets](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/RFBucket_Matching.ipynb)
* [deformation of RF bucket potential by longitudinal space charge](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/RFBucket_Potential_Deformation.ipynb)
* [radio frequency quadrupoles for Landau damping](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/RFQTest.ipynb)
* [slicing in the longitudinal plane](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/SlicingTest.ipynb)
* [space charge models on the GPU](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/Space_Charge_on_the_GPU.ipynb)
* [Bassetti-Erskine / Gaussian space charge model](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/TransverseGaussianSpaceCharge.ipynb)
* [transverse tracking](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/TransverseTrackingTest.ipynb)
* [wakefield models with dipolar and quadrupolar terms](https://nbviewer.jupyter.org/github/PyCOMPLETE/PyHEADTAIL-playground/blob/master/howto_notebooks/WakeTest.ipynb)

# Further Resources:
The [PyHEADTAIL wiki](https://github.com/PyCOMPLETE/PyHEADTAIL/wiki) is the entry level for documentation on PyHEADTAIL.

These are some potentially useful and concise slides as an overview on PyHEADTAIL, some from the [PyHEADTAIL developer meetings](https://indico.cern.ch/category/6360/):
* [PyHEADTAIL Playground, testing and release process](https://indico.cern.ch/event/799991/)
* [GPU integration in PyHEADTAIL](https://indico.cern.ch/event/807697/)
* [Space Charge integration in PyHEADTAIL](https://indico.cern.ch/event/828068/)
