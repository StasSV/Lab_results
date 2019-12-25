import numpy as np
from decimal import *
from random import random, choice
import matplotlib.pyplot as plt
import os.path
import math
import xlwt

MAX_iterations = 4000
population = []
health = []
NI = 0
length = 0
pop_number = 0


def progon_rws(pop_l, pop_num, coef_increase, increse_iter, pm, progon_number):
    global population
    global NI
    global length
    global pop_number

    population = []
    NI = 0
    length = 0
    pop_number = 0

    pop = createPop(pop_l, pop_num)
    polimorf = 0
    polimorf_iter = 0
    check_polimorf = True
    excel_count = 1

    wb = xlwt.Workbook()
    sheet_hamming = wb.add_sheet('Hamming')
    sheet_pair = wb.add_sheet('Pair')
    sheet_wild = wb.add_sheet('Wild')

    for i in range(MAX_iterations):
        still_increase = i < increse_iter
        [chosen_pop, chosen_health] = rws(pop, coef_increase, still_increase)
        px = 1 / (10 * length)
        mutated_pop = mutation(chosen_pop, px, pm)
        print(i)
        if (i % 20 == 0) or (still_increase is True):
            print('draw histos')
            polimorf = countPolimorf(mutated_pop)
            if (polimorf >= 0.38) and (check_polimorf is True):
                polimorf_iter = i
                check_polimorf = False
            hd = hamming_distance(mutated_pop)
            pd = pair_distance(mutated_pop)
            wd = dist_to_wild(mutated_pop)
            hd_values = list(hd.values())
            pd_values = list(pd.values())
            wd_values = list(wd.values())
            sheet_hamming.write(0, excel_count, i)
            sheet_pair.write(0, excel_count, i)
            sheet_wild.write(0, excel_count, i)
            for k in range(1, length + 1):
                sheet_hamming.write(k, excel_count, hd_values[k - 1])
                sheet_pair.write(k, excel_count, pd_values[k - 1])
                sheet_wild.write(k, excel_count, wd_values[k - 1])
            excel_count += 1
            build_histogram(i, hd, 'rws', coef_increase, increse_iter, pop_num, 'hamming', progon_number, polimorf, pm)
            build_histogram(i, pd, 'rws', coef_increase, increse_iter, pop_num, 'pair', progon_number, polimorf, pm)
            build_histogram(i, wd, 'rws', coef_increase, increse_iter, pop_num, 'wild', progon_number, polimorf, pm)

        pop = mutated_pop

        if i == MAX_iterations - 1:
            print('stop')
            NI = i
            polimorf = countPolimorf(mutated_pop)
            break

    f = open("Task2 l={0}_N={1}_pm={2}_growthCoef={3}_growthIter={4}_selection={5}_P={6}.txt"
             .format(pop_l, pop_num, pm, coef_increase, increse_iter, 'rws', progon_number), "a")
    to_write = "1: {0} 2: {1} 3: {2}".format(NI, polimorf, polimorf_iter)
    f.write(to_write)
    f.close()
    wb.save("l={0}_N={1}_pm={2}_growthCoef={3}_growthIter={4}_selection={5}_P={6}.xls"
            .format(pop_l, pop_num, pm, coef_increase, increse_iter, 'rws health l', progon_number))
    return


def progon_tournament(pop_l, pop_num, t, coef_increase, increse_iter, pm, progon_number):
    global population
    global NI
    global length
    global pop_number

    population = []
    NI = 0
    length = 0
    pop_number = 0

    pop = createPop(pop_l, pop_num)
    polimorf = 0
    polimorf_iter = 0

    tour = 'tournament ' + str(t)

    for i in range(MAX_iterations):
        still_increase = i < increse_iter
        chosen_pop = tournament(pop, t, coef_increase, still_increase)
        px = 1 / (10 * length)
        mutated_pop = mutation(chosen_pop, px, pm)
        if i % 10 == 0:
            print('draw histos')
            polimorf = countPolimorf(mutated_pop)
            if polimorf >= 0.38:
                polimorf_iter = i
            hd = hamming_distance(mutated_pop)
            pd = pair_distance(mutated_pop)
            wd = dist_to_wild(mutated_pop)
            build_histogram(i, hd, tour, coef_increase, increse_iter, pop_num, 'hamming', progon_number, polimorf)
            build_histogram(i, pd, tour, coef_increase, increse_iter, pop_num, 'pair', progon_number, polimorf)
            build_histogram(i, wd, tour, coef_increase, increse_iter, pop_num, 'wild', progon_number, polimorf)

        pop = mutated_pop

        if i == MAX_iterations - 1:
            NI = i
            polimorf = countPolimorf(mutated_pop)
            break

    f = open("Task2 l={0}_N={1}_pm={2}_growthCoef={3}_growthIter={4}_selection={5}_P={6}.txt"
             .format(pop_l, pop_num, pm, coef_increase, increse_iter, tour, progon_number), "a")
    to_write = "1: {0} 2: {1} 3: {2}".format(NI, polimorf, polimorf_iter)
    f.write(to_write)
    f.close()
    return


# create population of N with length l
def createPop(l, pop_N):
    global population
    global length
    global pop_number
    length = l
    pop_number = pop_N
    for i in range(pop_N):
        population.append([])
        for j in range(l):
            population[i].append(0)

    return population


# return array of health
def calcHealth(pop):
    for i, val in enumerate(pop):
        health.append(calchGenHeath(val))
    return health


def calchGenHeath(gen):
    return len(gen)


# define pm and mutate
def mutation(pop, px, function):
    return {
        function == 1: mutate(pop, px),
        function == 2: mutate(pop, px + 0.2 * px),
        function == 3: mutate(pop, px - 0.2 * px),
        function == 4: mutate(pop, px / 2),
        function == 5: mutate(pop, px / 10),
        function == 6: mutate(pop, px / 100)
    }[True]


def mutate(pop, pm):
    for i, val in enumerate(pop):
        for j, loc in enumerate(val):
            u = random()
            if u >= pm:
                if loc == 0:
                    val[j] = 1
                else:
                    val[j] = 0
    return pop


def rws(pop, coef_increase, still_increase):
    popHealth = calcHealth(pop)
    sum_health = sum(popHealth)
    pop = sorted(pop, key=calchGenHeath)
    probability = []
    chosen = []
    chosen_health = []
    for i in pop:
        # prob = Decimal(calchGenHeath(i))/Decimal(sum_health)
        prob = calchGenHeath(i) / sum_health
        probability.append(prob)
    cumsum = np.cumsum(probability)
    cumsum[len(cumsum) - 1] = Decimal(1.)

    if not still_increase:
        coef_increase = 1

    for i in range(math.ceil(coef_increase * len(pop))):
        u = random()
        for j, val in enumerate(cumsum):
            if val >= u:
                chosen_health.append(calchGenHeath(pop[j]))
                chosen.append(pop[j].copy())
                break
    global pop_number
    pop_number = len(chosen)
    return chosen, chosen_health


def selRandom(individuals, k):
    return [choice(individuals) for i in range(k)]


def tournament(pop, k, coef_increase, still_increase):
    chosen = []
    if not still_increase:
        coef_increase = 1

    for i in range(math.ceil(len(pop) * coef_increase)):
        aspirants = selRandom(pop, k)  # select k and select the best of them
        chosen.append(max(aspirants, key=calchGenHeath))

    global pop_number
    pop_number = len(chosen)
    return chosen


def countPolimorf(pop):
    polimorf = 0
    for i in range(length):
        for j in range(pop_number):
            if pop[j][i] == 1:
                polimorf += 1
                break
    return polimorf / length


def hamming_distance(pop):
    hammings = {}
    for i in range(length + 1):
        hammings.update({'%i' % i: 0})

    for i in range(len(pop)):
        count_diff = 0
        for j, gen in enumerate(pop[i]):
            if gen == 1:
                count_diff += 1

        if '%d' % count_diff in hammings:
            hammings['%d' % count_diff] += 1

    return hammings


def dist_to_wild(pop):
    wild_gen = []
    distances = {}
    for i in range(length + 1):
        distances.update({'%i' % i: 0})

    for i in range(length):
        ones = 0
        zeros = 0
        for j in range(pop_number):
            if pop[j][i] == 0:
                ones += 1
            elif pop[j][i] == 1:
                zeros += 1
        if ones > zeros:
            wild_gen.append(0)
        else:
            wild_gen.append(1)

    for i in range(len(pop)):
        count_diff = 0
        for j, gen in enumerate(pop[i]):
            if gen != wild_gen[j]:
                count_diff += 1

        if '%d' % count_diff in distances:
            distances['%d' % count_diff] += 1

    return distances


def pair_distance(pop):
    distances = {}
    for i in range(length + 1):
        distances.update({'%i' % i: 0})

    for i, gen in enumerate(pop):
        for j in range(i + 1, len(pop)):
            comp_gen = pop[j]
            count_diff = 0
            for k in range(length):
                if gen[k] != comp_gen[k]:
                    count_diff += 1
            if '%d' % count_diff in distances:
                distances['%d' % count_diff] += 1

    return distances


# Save a histogram to png file with all the parameters specified (except mutation)
def build_histogram(iter_num, distances, selection_type, coef_increase, iter_increase, pop_num, type, progon, polimorf,
                    pm):
    current_directory = os.getcwd() + '/Plots2'
    results_dir = os.path.join(current_directory,
                               'Selection={0},CoefIncrease={1},IterIncrese={2},N={3};l={4};type={5};pm={6};Progon={7}'
                               .format(selection_type, coef_increase, iter_increase, pop_num, length, type, pm, progon))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    plt.title("Polimorf: {0}%".format(polimorf * 100))
    plt.bar(list(distances.keys()), distances.values(), color='g', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig(results_dir + "/N={0}_l={1}_iter={2}.png".format(pop_number, length, iter_num))


# progon_rws(pop_l, pop_num, coef_increase, increse_iter, pm, progon_number)
progon_rws(10, 1, 2, 10, 1, 1)
progon_rws(10, 1, 2, 10, 2, 1)
progon_rws(10, 1, 2, 10, 3, 1)
progon_rws(10, 1, 2, 10, 4, 1)
progon_rws(10, 1, 2, 10, 5, 1)
progon_rws(10, 1, 2, 10, 6, 1)

progon_rws(10, 1, 1.005, 3500, 1, 1)
progon_rws(10, 1, 1.005, 3500, 2, 1)
progon_rws(10, 1, 1.005, 3500, 3, 1)
progon_rws(10, 1, 1.005, 3500, 4, 1)
progon_rws(10, 1, 1.005, 3500, 5, 1)
progon_rws(10, 1, 1.005, 3500, 6, 1)

progon_rws(10, 10, 2, 10, 1, 1)
progon_rws(10, 10, 2, 10, 2, 1)
progon_rws(10, 10, 2, 10, 3, 1)
progon_rws(10, 10, 2, 10, 4, 1)
progon_rws(10, 10, 2, 10, 5, 1)
progon_rws(10, 10, 2, 10, 6, 1)

progon_rws(10, 10, 1.005, 3500, 1, 1)
progon_rws(10, 10, 1.005, 3500, 2, 1)
progon_rws(10, 10, 1.005, 3500, 3, 1)
progon_rws(10, 10, 1.005, 3500, 4, 1)
progon_rws(10, 10, 1.005, 3500, 5, 1)
progon_rws(10, 10, 1.005, 3500, 6, 1)
