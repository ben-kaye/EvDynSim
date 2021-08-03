import math
# FUNCTIONS FOR EV SIM


def select_bact(pops, DFE, dt, cull_thresh):
    # implement formula for fixed population

    new_pops = [get_growth(p, dt, DFE) for p in pops]
    pop_total = sum(new_pops)

    for i in range(len(new_pops)):
        f_new = new_pops[i]/pop_total

        if f_new < cull_thresh:
            pops[i][1] = 0
        else:
            pops[i][1] = f_new  # assign normalised pop


def mut_bact(pops, DFE, PDF, tot_mut_rate, RNG):
    
    ## KEY parameter ##
    exp_size = 30

    for i in range(len(pops)):
        bpop = pops[i]
        exp_m = bpop[1]*tot_mut_rate  # get expected value
        mut_n = int(RNG.poisson(exp_m))

        for k in range(mut_n):
            new_genome = list(bpop[0])  # break reference
            new_genome.append(get_mut(PDF, DFE, RNG))

            new_genome.sort()

            
            size = RNG.poisson(exp_size)

            freq = size*1e-6  # pretend total pop is at 1,000,000

            bpop[1] -= freq

            # not sure if want to cast to list
            existing_gens = list(map(lambda u: u[0], pops))
            if new_genome in existing_gens:
                index = existing_gens.index(new_genome)
                pops[index][1] += freq
            else:
                pops.append([new_genome, freq])


def get_PDF(DFE):

    tot_rate = mut_rate_total(DFE)

    PDF = list(map(lambda x: x[2]/tot_rate, DFE.values()))

    return PDF


def mut_rate_total(DFE):
    return sum(map(lambda p: p[2], DFE.values()))


def get_mut(PDF, DFE, RNG):
    return RNG.choice(list(DFE), p=PDF, replace=False)

def get_fitness(pop, DFE):
    return sum(map(lambda g: DFE[g][1], pop[0]))

def augmented_fitness(pop, DFE):
    # define m
    m = 0.3

    # return s_aug = s_0 - m*susceptiblity
    return sum(map(lambda g: DFE[g][1] - m*DFE[g][3], pop[0]))

def selection_prob(chi):
    p = chi/0.6
    if p > 1:
        p = 1
    return p


def get_growth(pop, dt, DFE):
#TODO FIX RETENTION CALC
    fitness_total = get_fitness(pop, DFE)
    retention_coeff = selection_prob(DFE[pop[0][-1]][2]) # TODO need to fix, currently takes last magnetism probs want biggest
    QbyV = 1e-2  # per gen

    return pop[1]*math.exp(dt*(fitness_total + retention_coeff*QbyV))
