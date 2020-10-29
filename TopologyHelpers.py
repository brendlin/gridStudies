# Helper functions to return topologically-relevant output
# Input includes grid2op environment

import CommonHelpers
import numpy as np

def MakeAdjacencyMatrix(_env,n_buses=2,skipExternals=False,printLineIDs=False,lineOnBits=-1) :
    # Given an environment, make the adjacency matrix
    # (assuming all lines are on, and all buses are fully connected).

    n_entries = len(_env.sub_info)
    if not skipExternals : 
        n_entries += _env.n_gen + _env.n_load

    if n_buses == 2 :
        n_entries += len(_env.sub_info)

    fill_value = -1 if printLineIDs else 0
    adj_matrix = np.full(shape=[n_entries,n_entries],fill_value=fill_value)

    for i in range(_env.n_line) :

        if (lineOnBits >= 0) and not ((0b1 << i) & lineOnBits) :
            #print('line is turned off:',i)
            continue

        val = i if printLineIDs else 1
        bus = 0
        lor = n_buses*(_env.line_or_to_subid[i]) + bus
        lex = n_buses*(_env.line_ex_to_subid[i]) + bus
        adj_matrix[lor][lex] = val
        adj_matrix[lex][lor] = val

    if skipExternals :
        return adj_matrix

    for i in range(_env.n_gen) :
        bus = 0
        offset = n_buses*len(_env.sub_info)
        gi = i + offset
        sub = n_buses*(_env.gen_to_subid[i]) + bus
        adj_matrix[gi][sub] = 1
        adj_matrix[sub][gi] = 1

    for i in range(_env.n_load) :
        bus = 0
        offset = n_buses*len(_env.sub_info) + _env.n_gen
        li = i + offset
        sub = n_buses*(_env.load_to_subid[i]) + bus
        adj_matrix[li][sub] = 1
        adj_matrix[sub][li] = 1

    return adj_matrix


class AdjacencyMatrixClass :

    def __init__(self,env) :
        self.adjacency_matrix = MakeAdjacencyMatrix(env,n_buses=2,skipExternals=False)
        self.n_sub = len(env.sub_info)
        self.n_gen = env.n_gen
        self.n_load = env.n_load

    def ElementIDToAdjacencyIndex(self,elementID) :
        isLineOr = (elementID//1000 == 1)
        isLineEx = (elementID//1000 == 2)
        isLoad = (elementID//1000 == 3)
        isGen = (elementID//1000 == 4)

        if isGen :
            return self.n_sub*2 + (elementID - 4000)
        if isLoad :
            return self.n_sub*2 + self.n_gen + (elementID - 3000)
        if isLineEx :
            return (elementID - 2000)*2
        return (elementID - 1000)*2


def MakeLaplacian(_env,n_buses=2,skipExternals=False,lineOnBits=-1) :
    # Given an environment, make the Laplacian matrix
    # (assuming all lines are on, and all buses are fully connected).

    lap_matrix = -1*MakeAdjacencyMatrix(_env,
                                        n_buses=n_buses,
                                        skipExternals=skipExternals,
                                        lineOnBits=lineOnBits)

    if skipExternals :
        for i,row in enumerate(lap_matrix) :
            lap_matrix[i][i] = -sum(row)
        return lap_matrix

    for i,info in enumerate(_env.sub_info) :
        lap_matrix[i*n_buses][i*n_buses] = info

    for i in range(_env.n_gen) :
        bus = 0
        offset = n_buses*len(_env.sub_info)
        lap_matrix[i+offset][i+offset] = 1

    for i in range(_env.n_load) :
        bus = 0
        offset = n_buses*len(_env.sub_info) + _env.n_gen
        lap_matrix[i+offset][i+offset] = 1

    return lap_matrix


def IsConnectedLaplacianEigenvalue(_lap_matrix) :
    evals,evecs = np.linalg.eig(_lap_matrix)
    #print('Eigenvalues/vectors:')
    #for i in evals :
    #    print('{:.01f} + {:.01f}j'.format(i.real,i.imag))

    n_zeroeig = 0
    for i in evals :
        if abs(i.real) < 1e-9 :
            n_zeroeig += 1
    #print('Number of zero eigenvalues:',n_zeroeig)

    return (n_zeroeig <= 1)

def IsConnectedManual(_adj_matrix) :
    i_start = 0
    SetOfAlreadyTraversed = set([i_start])

    new_vertices = [i_start]
    while True :
        next_new_vertices = []
        for i_vert in new_vertices :
            for j_vert,entry in enumerate(_adj_matrix[i_vert]) :
                if (entry == 1) and (j_vert not in SetOfAlreadyTraversed) :
                    SetOfAlreadyTraversed.add(j_vert)
                    next_new_vertices.append(j_vert)

        # check if there are no new vertices
        if not len(next_new_vertices) :
            break

        new_vertices = next_new_vertices

    return len(SetOfAlreadyTraversed) == len(_adj_matrix)


def GetDisjointSets(_adj_matrix) :
    # Return a dictionary of disjoint sets.

    i_start = 0

    # Consider the case where you might have >1 disjoint sets
    AllTraversedVertices = set([i_start])

    DisjointSetsOfAlreadyTraversed = dict()
    DisjointSetsOfAlreadyTraversed[i_start] = set([i_start])
    CurrentSetOfAlreadyTraversed = DisjointSetsOfAlreadyTraversed[i_start]

    new_vertices = [i_start]
    while True :

        #print('new_vertices:',new_vertices)
        #print('DisjointSetsOfAlreadyTraversed',DisjointSetsOfAlreadyTraversed)
        #print('All traversed vertices',AllTraversedVertices)

        next_new_vertices = []
        for i_vert in new_vertices :
            for j_vert,entry in enumerate(_adj_matrix[i_vert]) :
                if (entry == 1) and (j_vert not in CurrentSetOfAlreadyTraversed) :
                    CurrentSetOfAlreadyTraversed.add(j_vert)
                    AllTraversedVertices.add(j_vert)
                    next_new_vertices.append(j_vert)

        # check if there are no new vertices
        if not len(next_new_vertices) :

            # If you traversed every vertex, then break
            if len(CurrentSetOfAlreadyTraversed) == len(_adj_matrix) :
                break

            if len(AllTraversedVertices) == len(_adj_matrix) :
                break

            # If not, start a new Set with the next untouched vertex
            for i_vert in range(len(_adj_matrix)) :
                if i_vert in AllTraversedVertices :
                    continue
                next_new_vertices.append(i_vert)
                DisjointSetsOfAlreadyTraversed[i_vert] = set([i_vert])
                CurrentSetOfAlreadyTraversed = DisjointSetsOfAlreadyTraversed[i_vert]
                AllTraversedVertices.add(i_vert)
                break

        new_vertices = next_new_vertices

    return DisjointSetsOfAlreadyTraversed


# A little bit more parseable version
def PrintAdjacencyMatrix(adj_matrix,skipExternals=False,nullstr='·') :
    tmp = '     '
    for i in range(len(env.sub_info)) :
        tmp += 's%02d '%(i)
    print(tmp)

    nDisplay = 2*len(env.sub_info) + env.n_gen + env.n_load
    if skipExternals :
        nDisplay = 2*len(env.sub_info)

    for index,i in enumerate(adj_matrix[:nDisplay]) :
        tmp = ''
        if index < 2*len(env.sub_info) :
            tmp += ('s%02d '%(index/2) if not index%2 else '    ')
        else :
            tmp += '    '

        for jndex, j in enumerate(i[:nDisplay]) :
            tmp_nullstr = nullstr
            if jndex >= 2*len(env.sub_info) and jndex < 2*len(env.sub_info) + env.n_gen :
                tmp_nullstr = '·'
            if index >= 2*len(env.sub_info) and index < 2*len(env.sub_info) + env.n_gen :
                tmp_nullstr = '·'
            tmp += ' ' + ('%d'%(j) if j>=0 else tmp_nullstr)
        print(tmp)

    return


# Print the line IDs
def PrintLineIDs(_env) :

    line_lookup = MakeAdjacencyMatrix(_env,n_buses=1,skipExternals=True,printLineIDs=True)

    tmp = '     '
    for i in range(len(_env.sub_info)) :
        tmp += 'sub%02d '%(i)
    print(tmp)

    for i,mline in enumerate(line_lookup[:2*len(_env.sub_info)]) :        
        tmp = ''
        tmp += ('sub%02d '%(i/2) if not i%2 else '      ')
        tmp += ' '.join(('%02d'%(a) if a>=0 else '· ') for a in list(mline[:2*len(_env.sub_info)]))
        print(tmp)

    return


def nDisconnected(flag,n_total) :
    # Quick function: return the number of disconnected lines, given a bitwise flag
    flag_str = "{:0{}b}".format(flag,n_total)
    n_zeros = flag_str.count('0')
    return n_zeros


def nConnected(flag,n_total) :
    # Quick function: return the number of connected lines, given a bitwise flag
    flag_str = "{:0{}b}".format(flag,n_total)
    n_ones = flag_str.count('1')
    return n_ones


def ExcludedByBitsetWithFewerDisconnections(line_bitset,excluded_bitsets) :
    # If there is a bitset in the list "excluded_bitsets"
    # whose list of off-lines is a subset of the list of off-lines of this line_bitset,
    # then return True.

    all_on = CommonHelpers.FullyConnectedBitset(20)
    line_bitset_flipped = all_on - line_bitset

    for excl_bitset in excluded_bitsets :

        excl_bitset_flipped = all_on - excl_bitset

        #print('Checking 0b{:b} against 0b{:b}'.format(line_bitset_flipped,excl_bitset_flipped))
        if (line_bitset_flipped & excl_bitset_flipped) == excl_bitset_flipped :
            #print('Checking 0b{:b} against 0b{:b} - excluded.'.format(line_bitset,excl_bitset))
            return True

    return False


def AddExcludedBitset(line_bitset,excluded_bitsets) :
    # If "excluded_bitsets" is a list of already-excluded (by topology) bitsets,
    # then add this excluded bitset to that list, with a caveat:
    #  - If there is another bitset, already in this list, that is excluded, where the off-lines
    #    of one bitset is a subset of the list of off-lines of the other bitset, then only keep the
    #    bitset with a smaller list of off-lines (since the other is definitely not a
    #    minimum-cut bitset).

    all_on = CommonHelpers.FullyConnectedBitset(20)
    line_bitset_flipped = all_on - line_bitset

    for i in range(len(excluded_bitsets)-1,-1,-1) :
        excl_bitset = excluded_bitsets[i]
        excl_bitset_flipped = all_on - excl_bitset

        if (line_bitset_flipped & excl_bitset_flipped) == excl_bitset_flipped :
            # Already covered - do nothing
            print('0b{:b} already covered by 0b{:b}'.format(excl_bitset_flipped,line_bitset_flipped))
            return

        if (line_bitset_flipped & excl_bitset_flipped) == line_bitset_flipped :
            # pop the other bitset, add this new bitset and exit
            #print('0b{:b} is superseded by 0b{:b}'.format(excl_bitset_flipped,line_bitset_flipped))
            excluded_bitsets.pop(i)

    # No conflict; add line_bitset
    excluded_bitsets.append(line_bitset)
    #print('Excluded bitsets changed; new size is',len(excluded_bitsets))

    return
