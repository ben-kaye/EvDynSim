import math
# FUNCTIONS FOR EV SIM


def select_bact(pops, DFE, dt, cull_thresh):
    # logistic growth + eliminate pop with size < 10
    s_bar = get_fitness_avg(pops, DFE)

    for bpop in pops:
        s_tot = get_fitness(bpop, DFE)
        if bpop[1] > 0:
            bpop[1] *= math.exp(dt*(s_tot - s_bar))

        bpop[1] = int(bpop[1])

        if bpop[1] <= cull_thresh:
            bpop[1] = 0


def mut_bact(pops, DFE, PDF, tot_mut_rate, RNG):
    for i in range(len(pops)):
        bpop = pops[i]
        exp_m = bpop[1]*tot_mut_rate  # get expected value
        mut_n = int(RNG.poisson(exp_m))

        for k in range(mut_n):
            new_genome = list(bpop[0])  # break reference
            new_genome.append(get_mut(PDF, DFE, RNG))

            new_genome.sort()

            exp_size = 30
            size = RNG.poisson(exp_size)

            bpop[1] -= size

            # not sure if want to cast to list
            existing_gens = list(map(lambda u: u[0], pops))
            if new_genome in existing_gens:
                index = existing_gens.index(new_genome)
                pops[index][1] += size
            else:
                pops.append([new_genome, size])


def get_PDF(DFE):

    tot_rate = mut_rate_total(DFE)

    PDF = list(map(lambda x: x[2]/tot_rate, DFE.values()))

    return PDF


def mut_rate_total(DFE):
    return sum(map(lambda p: p[2], DFE.values()))


def get_mut(PDF, DFE, RNG):
    return RNG.choice(list(DFE), p=PDF, replace=False)


def get_fitness_avg(pops, DFE):
    total = sum(map(lambda u: u[1], pops))
    # returns list of fits
    fits = map(lambda h: get_fitness(h, DFE), pops)

    # return avg fitness
    return sum(map(lambda u, v: u[1]*v, pops, fits))/total


def get_fitness(pop, DFE):
    return sum(map(lambda g: DFE[g][1], pop[0]))
