import fixedN_ev as EV
from numpy import random as RND


import time

# data structures:
# DFE // id: [name, metabolic_advantage, mut_rate, susceptibility]
# pops // [[[mutations], pop_size] ... ]


### PARAMS ###

n_gens = int(5e5) 
cull_thresh = 1e-6 # 1 cell per million
dt = 1 # in units generation

Npop = 1e7 # affects mutation count (NU)

nNratio = 1e-2 # just for binary pop test


# thought I might bundle up the parameters, but no good data structure for this
#params = {}

DFE = {
    0: ['WT', 0, 0, 0],
    1: ['fur', -2e-3, 1e-7, 1e-2],
    2: ['fur2', -1e-3, 5e-7, 4e-2]
}

###      ###

PDF = EV.get_PDF(DFE)
tot_mut_rate = EV.mut_rate_total(DFE)
RNG = RND.default_rng()

pops = [[[0], (1 - nNratio)], [[0,1], nNratio]]

alt_pops = list(pops)

start_time = time.time()

### SIMULATION LOOP ###
for k in range(1, n_gens):

    # grow population and take away
    EV.select_bact(pops, DFE, dt, cull_thresh)

    # mutate
    EV.mut_bact(pops, DFE, PDF, Npop, tot_mut_rate, RNG)

elapsed_time = time.time() - start_time

print('program finished, exec time: ', "{:.3f}".format(elapsed_time), 's')
