# DynGillEpi

Python wrapper for the Gillespie contagion functions developed and written by CL Vestergaard and M Génois in [Temporal Gillespie Algorithm: Fast Simulation of Contagion Processes on Time-Varying Networks](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004579).
Check out the original source code at https://github.com/CLVestergaard/TemporalGillespieAlgorithm.

Note that I changed the source code in the way that all random number generators are now taken from the C++ standard library.

## Install

### Python

    $ sudo pip install ./DynGillEpi

## Examples

### sandbox

    $ cd sandbox
    $ python SIS.py

### SIS.py

```python
import numpy as np
import matplotlib.pyplot as pl
import networkx as nx
import DynGillEpi

SIS = DynGillEpi.SIS_Poisson_homogeneous

N_nodes = 10
n_slices = 10
contact_list = [ [tuple(e) for e in nx.fast_gnp_random_graph(N_nodes,2./(N_nodes-1.)).edges()] for i in range(n_slices) ]


n_simulations = 100
T_simulation = 100

infection_rate = 10. # number of events per SI-link per dt
recovery_rate = 1.   # number of recoveries per infected per dt


result = SIS(N_nodes,
             contact_list,
             infection_rate,
             recovery_rate,
             T_simulation,
             number_of_simulations = n_simulations,
             initial_number_of_infected = 3,
             seed = 324345,
             verbose = False,
             )

all_I = np.array(result.I,dtype=float)
all_SI = np.array(result.SI,dtype=float)

obs = [all_I, all_SI]
t = np.arange(T_simulation)

for o in obs:
    mean = o.mean(axis=0)
    std = o.std(axis=0) / np.sqrt(n_simulations - 1.)
    pl.errorbar(t,mean,std)

pl.xlabel(r'time $t/\Delta t$')
pl.ylabel(r'total number')
pl.legend(['of infected', 'of SI-links'])
pl.show()

```
