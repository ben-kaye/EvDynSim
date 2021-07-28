import ev2
from numpy import random as RND


import time

# data structures:
# DFE // id: [name, selective_advantage, mut_rate]
# pops // [[[mutations], pop_size] ... ]


### PARAMS ###

start_pop_size = int(1e7)
n_gens = int(1e3)
cull_thresh = 10
dt = 1

DFE = {
    0: ['WT', 0, 0],
    1: ['PDE2', 8e-2, 2e-7],
    2: ['IRA1', 10e-2, 0.8e-7]
}

###      ###

PDF = ev2.get_PDF(DFE)
tot_mut_rate = ev2.mut_rate_total(DFE)
RNG = RND.default_rng()

pops = [[[0], start_pop_size]]

start_time = time.time()

for k in range(1, n_gens):
    ev2.mut_bact(pops, DFE, PDF, tot_mut_rate, RNG)
    ev2.select_bact(pops, DFE, dt, cull_thresh)

elapsed_time = time.time() - start_time

print('program finished, exec time: ', "{:.3f}".format(elapsed_time), 's')
