# mEC Model
 
 This code was written using NetPyNE (Networks with Python and NEURON). For more information about the tool, see [NetPyNE docs](http://doc.netpyne.org/)

## Download the Repository

To clone the repository, open a terminal in the directory where you'd like to store the project and run:

```bash
git clone https://github.com/RomanB22/mEC_model2025.git
````

After that, make sure to move to the repo folder using `cd mEC_model2025`.

## Installation

Install the required Python packages using:

```bash
pip install -r requirements.txt
```

## Simulations Setup

Add the root directory to `PYTHONPATH`, and compile the mechanisms:

```bash
export PYTHONPATH=PYTHONPATH:$PWD
nrnivmodl mod
```

### Setting up the parameters

In `src/cfg.py` you can modify all the parameters, variables and conditions for the simulation. For example, you can change the values of the flags to have sinusoidal optogenetic drive, choose the stellate cell model, add gap junction connectivity between the PV cells, change the recorded cells, simulation time, etc.

## ▶Running Simulations

### Run a Single Simulation

All simulation parameters are defined in the `src/cfg.py` file. To run a single simulation:

```bash
python -u src/init.py
```

### Run Batch Simulations (Experimental)

To run multiple simulations on the terminal using the batch submitter:

```bash
python -u src/batch.py
```
The function `runNetworks(NumNetworks=X)` will run X number of networks with different connectivity but the same parameters as in `src/cfg.py`. It is possible to create another function with another set of 
```bash
params['property'] = [list of values]
```
where that property is any of the ones defined in `src/cfg.py`. That will create automatically all possible combination of values between all the params (be careful since the number of combinations grows exponentially). See for example the function `runDifferentReversalPot()`

> ⚠️ Note: Batch simulations have not been fully tested yet.

## Model Description

*(Add your model description here.)*

## To-Do List

- [ ] Check that the function `SC_Mittal(cwd, cfg)` in `src/defs.py` is loading the correct cell (should be the same number as in the end of `Noise_Osc_july.hoc`)
- [ ] Once that is checked, we can modify the parameters of the SC for each model in the python code `src/defs.py`, and we can load all the different models from the .csv file
- [ ] The next step will be to modify the code to add heterogeneity in the SC population
- [ ] Check that the temperature for the simulation is correct (`cfg.hParams = {'celsius': 23, 'v_init': -80}` in `src/cfg.py`). That shouldn't change anything unless any of the ion channels have a temperature-dependent behavior
- [ ] Add instructions on how to run NEURON with `mpi` to speed up sims. Ask Chris since it depends on each cluster configuration.
- [x] Completed Item