#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
#
# linear optimisation to calculate minimum and maximum commute
#
# (C) Timotheus Andreas Klein 2020 - Tim.Klein@gmx.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import os
import numpy
import glob
from ipfn import ipfn
from pulp import *		# requires https://github.com/coin-or/pulp with mip start


# folders and
datapath = 'C:\\Users\\trist\\Documents\\veroeffentlichung\\Pendler\\daten\\lin_opt\\'
scriptpath = 'C:\\Program Files (x86)\\Python37-32\\Lib\\site-packages\\tims_sandbox\\'
        
# functions
# requires csv table format without header

def commuter_evaluation(ags):
    log('\nstarting evaluation for ags ' + ags)
    # sources
    # 1. costfile: residence key;employment key;direct distance;routed distance;commutes
    costfile = datapath + 'pg_out_regmatrix_'+ags+'.csv'
    # 2. row sums: total of residents in area
    residentfile = datapath + 'pg_out_regmatrix_'+ags+'_o.csv'
    # 3. column sums: total of workplaces in area
    workplacefile = datapath + 'pg_out_regmatrix_'+ags+'_d.csv'
    # compile input data
    # costs
    cost_analysis = read_matrix_lbl(costfile, 3)
    costs = cost_analysis[0]
    # commuter flows
    commute_analysis = read_matrix_lbl(costfile, 4)
    commute = commute_analysis[0]
    # zones
    zones = sorted(cost_analysis[1].keys())
    zone_index = indexzones(zones)
    # put costs into numpy matrix
    costmatrix = create_matrix_from_dict(costs, zones, zone_index[0])
    np_costs = numpy.array(costmatrix)
    # put commuters into numpy matrix
    commutematrix = create_matrix_from_dict(commute, zones, zone_index[0])
    np_commute = numpy.array(commutematrix)
    # read row- and columnsums: residents and workplaces
    residents = get_cmatrix_sums(residentfile, zones)
    workplaces = get_cmatrix_sums(workplacefile, zones)
    log('Sum residents:\t' + str(sum(residents.values())))
    log('Sum workplaces:\t' + str(sum(workplaces.values())))
    log('Matrix size: ' + str(numpy.shape(np_commute)))
    # output
    results = []
    # observed commute
    c_obs = sum(sum(np_costs*np_commute))
    log('Observed Commute:\t\t\t' + str(c_obs))
    log(' Commutes:\t' + str(sum(sum(np_commute))))
    log(' non-zeros in commute matrix:\t' + str(numpy.count_nonzero(np_commute)))
    results.append(c_obs)   # results[0]: c_obs
    # minimum commute: np_commute, np_costs, residents, workplaces, zone_index
    c_min = calculate_minimum_commute(np_commute, np_costs, residents, workplaces, zone_index)
    # control and export results
    if c_min[4] == 1:
        # control
        print('resident sums, if diverging:')
        examples = 0
        for z in zones:
            if c_min[1][z] - residents[z] != 0 and examples < 11:
                print(z, ': ', c_min[1][z], residents[z])
                examples += 1
        print('workplace sums, if diverging:')
        examples = 0
        for z in zones:
            if c_min[2][z] - workplaces[z] != 0 and examples < 11:
                print(z, ': ', c_min[2][z], workplaces[z])
                examples += 1
        # export
        resultmatrix = create_export_array(c_min[0], zone_index[1])
        csvexport = array2csv(resultmatrix, ';')
        resultfile = datapath + 'c_min_'+ags+'.csv'
        saveFile(resultfile, csvexport)
        results.append(c_min[3])    # results[1]: c_min
        print('Minimum Commute:\t\t'+str(sum(sum(np_costs*c_min[0]))))
        log(' Commutes:\t' + str(sum(sum(c_min[0]))))
        log(' non-zeros in min commute matrix:\t' + str(numpy.count_nonzero(c_min[0])))
    else:
        results.append(0)
    # maximum commute: np_commute, np_costs, residents, workplaces, zone_index
    c_max = calculate_maximum_commute(np_commute, np_costs, residents, workplaces, zone_index)
    # control and export results
    if c_max[4] == 1:
        # control
        print('resident sums, if diverging:')
        examples = 0
        for z in zones:
            if c_max[1][z] - residents[z] != 0 and examples < 11:
                print(z, ': ', c_max[1][z], residents[z])
                examples += 1
        print('workplace sums, if diverging:')
        examples = 0
        for z in zones:
            if c_max[2][z] - workplaces[z] != 0 and examples < 11:
                print(z, ': ', c_max[2][z], workplaces[z])
                examples += 1
        # export
        resultmatrix = create_export_array(c_max[0], zone_index[1])
        csvexport = array2csv(resultmatrix, ';')
        resultfile = datapath + 'c_max_'+ags+'.csv'
        saveFile(resultfile, csvexport)
        results.append(c_max[3])    # results[2]: c_max
        print('Maximum Commute:\t\t'+str(sum(sum(np_costs*c_max[0]))))
        log(' Commutes:\t' + str(sum(sum(c_max[0]))))
        log(' non-zeros in max commute matrix:\t' + str(numpy.count_nonzero(c_max[0])))
    else:
        results.append(0)
    # random commute: conversion, residents, workplaces, zone_index
    c_rnd = calculate_random_commute(0.00001, np_costs, residents, workplaces, zone_index, ags)
    log(str(c_rnd[2]) + ' iterations')
    log('Random Commute:\t\t\t\t' + str(c_rnd[1]))
    log(' avg non-zeros in rnd commute matrix:\t' + str(c_rnd[0]))
    results.append(c_rnd[1])    # results[3]: c_rnd
    return results

def calculate_minimum_commute(np_commute, np_costs, residents, workplaces, zone_index):
    # np_commute:   observed commutes
    # np_costs:     numpy array of costs
    # residents:    dictionary[zone]:=resident workers
    # workplaces:   dictionary[zone]:=workplaces
    # zone_index:   [0] zone2index, [1] index2zone dictionaries
    log('Minimum Commute:')
    zones = sorted(zone_index[0].keys())
    zone2index = zone_index[0]
    index2zone = zone_index[1]
    # create problem variables
    LpProb = LpProblem("Minimum_Commute_Problem",LpMinimize)
    # define 'index sets'
    # matrix size:
    matrix_size = len(zones)
    I = range(matrix_size)
    J = range(matrix_size)
    # define parameter cost[i,j]
    cost = {}
    for i in I:
        for j in J:
            cost[i,j] = np_costs[i][j]
    # create problem variables: commuter volumes for all relations
    # m_commute[i,j]
    m_commute = LpVariable.dicts(name="m_commute", indexs=(I, J), lowBound=0, cat=LpInteger)
    # create objective function
    # objective
    LpProb += lpSum(cost[i,j] * m_commute[i][j] for i in I for j in J)
    # constraints
    # define residents[i], workplaces[j]
    constraint_residents = LpVariable.dicts(name="c_residents", indexs=(I))
    constraint_workplaces = LpVariable.dicts(name="c_workplaces", indexs=(J))
    # constraints residents, workplaces
    for i in I:
            constraint_residents[i] = residents[index2zone[i]]
            LpProb += constraint_residents[i] == lpSum(m_commute[i][j] for j in J), ""
    for j in J:
            constraint_workplaces[j] = workplaces[index2zone[j]]
            LpProb += constraint_workplaces[j] == lpSum(m_commute[i][j] for i in I), ""
    # solving process
    status = -3
    # without initial values
##    (status, LpProb) = start_solving(LpProb)
##    print('status: ', status)
    # set initial values
    if status != 1:
        log('setting initial values (observed commute)')
        for i in I:
            for j in J:
                m_commute[i][j].setInitialValue(np_commute[i][j])
        # solve with initial values
        (status, LpProb) = start_solving(LpProb)
        print('status: ' + str(LpStatus[status]))
    # preparing the output
    resultmatrix = numpy.zeros(shape=(matrix_size,matrix_size))
    residents2 = {}
    workplaces2 = {}
    resultvalue = value(LpProb.objective)
    try:
        for i in I:
            for j in J:
                if value(m_commute[i][j]) != 0:
                    resultmatrix[i][j] = value(m_commute[i][j])
        for i in I:
            residents2[index2zone[i]] = sum(resultmatrix[i])
        for j in J:
            workplaces2[index2zone[j]] = sum(resultmatrix[:])[j]
        log('status: ' + str(LpStatus[status]))   # The solution status 
        log('result minimisation:\t\t' + str(resultvalue))
    except:
        log('status: ' + str(LpStatus[status]))
    return resultmatrix, residents2, workplaces2, resultvalue, status

def calculate_maximum_commute(np_commute, np_costs, residents, workplaces, zone_index):
    # np_commute:   observed commutes
    # np_costs:     numpy array of costs
    # residents:    dictionary[zone]:=resident workers
    # workplaces:   dictionary[zone]:=workplaces
    # zone_index:   [0] zone2index, [1] index2zone dictionaries
    log('Maximum Commute:')
    zones = sorted(zone_index[0].keys())
    zone2index = zone_index[0]
    index2zone = zone_index[1]
    # create problem variables
    LpProb = LpProblem("Maximum_Commute_Problem",LpMaximize)
    # define 'index sets'
    # matrix size:
    matrix_size = len(zones)
    I = range(matrix_size)
    J = range(matrix_size)
    # define parameter cost[i,j]
    cost = {}
    for i in I:
        for j in J:
            cost[i,j] = np_costs[i][j]
    # create problem variables: commuter volumes for all relations
    # m_commute[i,j]
    m_commute = LpVariable.dicts(name="m_commute", indexs=(I, J), lowBound=0, cat=LpInteger)
    # create objective function
    # objective
    LpProb += lpSum(cost[i,j] * m_commute[i][j] for i in I for j in J)
    # constraints
    # define residents[i], workplaces[j]
    constraint_residents = LpVariable.dicts(name="c_residents", indexs=(I))
    constraint_workplaces = LpVariable.dicts(name="c_workplaces", indexs=(J))
    # constraints residents, workplaces
    for i in I:
            constraint_residents[i] = residents[index2zone[i]]
            LpProb += constraint_residents[i] == lpSum(m_commute[i][j] for j in J), ""
    for j in J:
            constraint_workplaces[j] = workplaces[index2zone[j]]
            LpProb += constraint_workplaces[j] == lpSum(m_commute[i][j] for i in I), ""
    # solving process
    status = -3
    # without initial values
##    (status, LpProb) = start_solving(LpProb)
##    print('status: ', status)
    # set initial values
    if status != 1:
        log('setting initial values (observed commute)')
        for i in I:
            for j in J:
                m_commute[i][j].setInitialValue(np_commute[i][j])
        # solve with initial values
        (status, LpProb) = start_solving(LpProb)
        print('status: ', status)
    # preparing the output
    resultmatrix = numpy.zeros(shape=(matrix_size,matrix_size))
    residents2 = {}
    workplaces2 = {}
    resultvalue = value(LpProb.objective)
    try:
        for i in I:
            for j in J:
                if value(m_commute[i][j]) != 0:
                    resultmatrix[i][j] = value(m_commute[i][j])
        for i in I:
            residents2[index2zone[i]] = sum(resultmatrix[i])
        for j in J:
            workplaces2[index2zone[j]] = sum(resultmatrix[:])[j]
        log('status: ' + str(LpStatus[status]))   # The solution status 
        log('result maximisation:\t\t' + str(resultvalue))
    except:
        log('status: ' + str(LpStatus[status]))
    return resultmatrix, residents2, workplaces2, resultvalue, status

def start_solving(LpProb):
    # solve the problem
    status = -3
    try:
        status = LpProb.solve(pulp.PULP_CBC_CMD())
        log('cbc worked')
    except:
        log('cbc did not work')
        try:
            status = LpProb.solve(GLPK_CMD(path = 'C:\\Program Files (x86)\\glpk-4.65\\w64\\glpsol.exe'))
            log('glpk worked')
        except:
            log('glpk did not work')
            try:
                status = LpProb.solve(pulp.CHOCO_CMD())
                log('choco worked')
            except:
                log('choco did not work')
    return status, LpProb

def calculate_random_commute(conversion, np_costs, residents, workplaces, zone_index, ags):
    # calculate random commute with proportional fitting of random commute table until
    # change of average of results goes below conversion value [0.010 = 1.0%]
    log('Random Commute:')
    zones = sorted(zone_index[0].keys())
    zone2index = zone_index[0]
    index2zone = zone_index[1]
    # matrix size:
    matrix_size = len(zones)
    # aggregates:
    resident_aggregates = numpy.zeros(matrix_size)
    workplace_aggregates = numpy.zeros(matrix_size)
    for i in range(matrix_size):
        resident_aggregates[i] = max(0.001, residents[index2zone[i]])
        workplace_aggregates[i] = max(0.001, workplaces[index2zone[i]])
    msum_target = round(sum(resident_aggregates))
    aggregates = [resident_aggregates, workplace_aggregates]
    # dimensions:
    dimensions=[[0],[1]]
    # calculation:
    # first iteration:
    steps = []
    non_zeros = []
    rnd_seed = numpy.random.randint(100, size=(matrix_size,matrix_size))
    IPF = ipfn.ipfn(rnd_seed, aggregates, dimensions, verbose=0, rate_tolerance=0.001)
    rnd_commute_matrix = IPF.iteration()
    rnd_commute_matrix = rnd_commute_matrix * msum_target/sum(sum(rnd_commute_matrix))
    steps.append(sum(sum(rnd_commute_matrix*np_costs)))
    non_zeros.append(numpy.count_nonzero(rnd_commute_matrix))
    # second iteration:
    rnd_seed = numpy.random.randint(100, size=(matrix_size,matrix_size))
    IPF = ipfn.ipfn(rnd_seed, aggregates, dimensions)
    rnd_commute_matrix = IPF.iteration()
    rnd_commute_matrix = rnd_commute_matrix * msum_target/sum(sum(rnd_commute_matrix))
    steps.append(sum(sum(rnd_commute_matrix*np_costs)))
    non_zeros.append(numpy.count_nonzero(rnd_commute_matrix))
    c_rnd_old = sum(steps[:-1])/len(steps[:-1])
    c_rnd_new = sum(steps)/len(steps)
    # loop:
    while abs(c_rnd_old-c_rnd_new)/c_rnd_old > conversion:
        rnd_seed = numpy.random.randint(100, size=(matrix_size,matrix_size))
        IPF = ipfn.ipfn(rnd_seed, aggregates, dimensions)
        rnd_commute_matrix = IPF.iteration()
##        # adjust result
##        rnd_commute_matrix = rnd_commute_matrix * msum_target/sum(sum(rnd_commute_matrix))
        steps.append(sum(sum(rnd_commute_matrix*np_costs)))
        non_zeros.append(numpy.count_nonzero(rnd_commute_matrix))
        c_rnd_old = sum(steps[:-1])/len(steps[:-1])
        c_rnd_new = sum(steps)/len(steps)
    log(' Commutes:\t' + str(sum(sum(rnd_commute_matrix))))
##    # write last result for control (can take days for matrices with size > 1000*1000)
##    log('Random Commute:\t\t\t' + str(c_rnd_new))
##    resultmatrix = create_export_array(rnd_commute_matrix, index2zone)
##    mcsv = array2csv(resultmatrix, ';')
##    mfilename = datapath + 'rnd_commute_' + ags + '.csv'
##    saveFile(mfilename, mcsv)
    return int(round(sum(non_zeros)/len(non_zeros),0)), c_rnd_new, len(steps)

# other auxiliary

def array2csv(A, delimiter):
    print('array2csv')
    csv = ''
    for line in A:
        for item in line:
            csv = csv + str(item) + delimiter
        csv = csv + '\n'
    return csv

def create_export_array(A, index2zone):
    # create list representation of array A - from, to, volume
    print('create_export_array')
    L = []
    I = range(len(A))
    J = range(len(A[:][0]))
    for i in I:
        for j in J:
            if A[i][j] != 0:
                line = [index2zone[i], index2zone[j], A[i][j]]
                L.append(line)
    return L

def create_matrix_from_dict(costs, zones, zone2index):
    # build a matrix from the relation-cost dictionary
    count = 0
    size = len(zones)
    matrix = numpy.zeros((size,size),float)
    for c in costs:
        i = c.split('_')[0]
        j = c.split('_')[1]
        matrix[zone2index[i]][zone2index[j]] = costs[c]
        count = count + costs[c]
    return matrix

def csv2array(csvdata, delimiter):
    A = []
    for line in csvdata.split('\n'):
        A.append(line.split(delimiter))
    return A

def get_cmatrix_sums(filename, allzones):
    # return commuter matrix sums from table prepared in postgres in dictionary:
    # zone;sum
    # no header
    sumsdict = {}
    for z in allzones:
        sumsdict[z] = 0
    sums_table = csv2array(readFile(filename), ';')
    for s in sums_table:
        try:
            sumsdict[s[0]] = float(s[1])
        except:
            pass
    return sumsdict

def indexzones(zones):
    # create index<>zone_number dictionaries
    zone2index = {}
    index2zone = {}
    for i, z in enumerate(zones):
        zone2index[z] = i
        index2zone[i] = z
    return zone2index, index2zone

def log(info):
    print(info)
    dirfile1str = datapath + 'log.txt'
    f = open(dirfile1str, 'a')
    f.write(info + '\n')
    f.close()

def read_matrix_lbl(matrixfile, index):
    # read matrix file line-by-line and store relational volumes in dictionary
    # customized for matrix file from;to;direct distance;routed distance;value
    # no header
    count = 0
    zones = {}
    matrixdict = {}
    f = open(matrixfile, 'r')
    for line in f:
        try:
            r = line.split(';')
            matrixdict[str(r[0])+'_'+str(r[1])] = float(r[index])
            zones[str(r[0])] = ''
            count = count + float(r[index])
        except:
            print(line)
            print(sys.exc_info())
            break
    log('read_matrix_lbl: matrix sum:\t' + str(count))
    return matrixdict, zones

def readFile(dirfile1str):
    f = open(dirfile1str, 'r')
    return f.read()

def saveFile(dirfile1str, info):
    f = open(dirfile1str, 'w')
    f.write(info)
    f.close()

# main
logfile = datapath + 'log.txt'
saveFile(logfile, 'Log:\n====')
result_dict = {}
ags_list = [
'12053000',
'12069397',
'14625020',
'12069656',
'11000000',
'12054000',
'12069616',
'12069590',
'12069304',
'12052000',
'12069604',
'14625040',
'12069454',
'12069596',
'12069017',
'12061448',
'14625250',
'05162022',
'14628270',
'15091375',
'03404000',
'14625480',
'05162016',
'05162020',
'15003000',
'05911000',
'14612000',
'14628060',
'15001000',
'05158028',
'14628110',
'05158012',
'14627210',
'14628400',
'05162012',
'05111000',
'14627010',
'14627060',
'05158032',
'03251037',
'03251048',
'03401000',
'14627140',
'03461006',
'03251047',
'05162024',
'03256014',
'05158036',
'15081045',
'03361013',
'14627230',
'05158024',
'04011000',
'13003000',
'05162008',
'05158004',
'03356007',
'15081105',
'03361009',
'13072057',
'13072098',
'03356003',
'05158008',
'03356002',
'03361008',
'03356011',
'14511000',
'05158016',
'05162004',
'14524330',
'05162028',
'05158026',
'15002000',
'14713000',
'09761000',
'14523320',
'03353018',
'16052000',
'01002000',
'09463000',
'09775135',
'08421000',
'16051000',
'01004000',
'16053000',
'07314000',
'07312000',
'06412000',
'07137060'
]

for ags in ags_list:
    result_dict[ags] = commuter_evaluation(ags)
    fileList = glob.glob(scriptpath+'*pulp.mps')
    for filename in fileList: os.remove(filename)
resultfilestr = ''
for k, v in result_dict.items():
    resultfilestr = resultfilestr + str(k) + ':' + str(v) + '\n'
resultfilename = datapath + 'results.csv'
saveFile(resultfilename, resultfilestr)
log(str(result_dict))

