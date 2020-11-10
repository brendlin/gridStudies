
import CommonHelpers
import BusTopologyHelpers
import numpy as np

class Substation :

    def __init__(self,id,nElements,elementIDs) :
        
        self.index = id
        self.nElements = nElements

        # Cache for keeping track of previously-checked bus states.
        # The key corresponds to the disconnected indices (binary number)
        self.validityCache = dict()

        # The nominal bus configuration is to put everything on bus 1
        self.currentBusConfig = CommonHelpers.FullyConnectedBitset(nElements)
        self.validityCache[self.currentBusConfig] = dict()

        # NEW type of element ID: 3 digits
        # First digit (starting from the left) is the type of element
        # Trailing three digits gives the position in the list (e.g. 0th line)
        # Examples:
        #  - 1000: Line ORIGIN, 0th in the list of lines
        #  - 2001: Line EXTREMILY, 1st in the list of lines
        #  - 3012: Load, 12th in the list of loads
        #  - 4005: Generator, 5th in the list of generators
        #
        # For example, the elements (in the order that they are described in the currentBusConfig)
        # could look like: [1000,1001,4004]
        # Meaning line origins #0 and #1 are connected, as is generator #4.
        # This should be static.
        self.elementIDs = np.array(elementIDs)

    # "Local" index here refers to the local bus bits here
    def LocalGeneratorIndices(self) :
        return list(np.where(self.elementIDs//1000 == 4)[0])
        
    def LocalLoadIndices(self) :
        return list(np.where(self.elementIDs//1000 == 3)[0])

    def LocalLineIndices(self) :
        return list(np.where(self.elementIDs//1000 <= 2)[0])
        
    def SetBusConfig(self,bits) :
        self.currentBusConfig = bits
        if bits not in self.validityCache.keys() :
            self.validityCache[bits] = dict()
        return

    def IsValidBooleanBusState(self,lineOnBits=-1) :

        localDisconnIndices = []
        cacheKey = 0

        if lineOnBits > 0 :
            for li in self.LocalLineIndices() :
                lineID = self.elementIDs[li]%1000
                isOn = lineOnBits & (0b1 << lineID)
                #print('{} is {}'.format(lineID,'on' if isOn else 'off'))
                if isOn :
                    continue
                else :
                    cacheKey += (0b1 << li)
                    localDisconnIndices.append(li)

        if cacheKey in self.validityCache[self.currentBusConfig].keys() :
            return self.validityCache[self.currentBusConfig][cacheKey]

        valid = BusTopologyHelpers.IsValidBooleanBusState(self.currentBusConfig,
                                                          self.nElements,
                                                          i_gens = self.LocalGeneratorIndices(),
                                                          i_loads = self.LocalLoadIndices(),
                                                          items_disconnected = localDisconnIndices,
                                                          verbose=False)

        self.validityCache[self.currentBusConfig][cacheKey] = valid
        return valid

    # def verifyBusConfiguration(self,adjacency_matrix_class) :
    #     for i in self.nElements :
    #         toAdjacencyID = list(adjacency_matrix_class.

def BuildSubstations(env) :

    # This builds substations from the environment, putting them into a format that is readily
    # useable by the tools developed to analyze available substation moves.

    sub_classes = []

    for sub in range(len(env.sub_info)) :

        # Element IDs, to be populated
        element_ids = [0]*env.sub_info[sub]

        # Find the line (OR) ids that link to this sub
        line_or_ids = np.where(env.line_or_to_subid == sub)[0]

        for lid in line_or_ids :

            # The sub position
            sub_pos = env.line_or_to_sub_pos[lid]
            if element_ids[sub_pos] != 0 :
                print('Error -- this sub position is already filled! (line OR)')

            element_ids[sub_pos] = 1000 + lid

        # Find the line (EX) ids that link to this sub
        line_ex_ids = np.where(env.line_ex_to_subid == sub)[0]
        for lid in line_ex_ids :
            sub_pos = env.line_ex_to_sub_pos[lid]
            if element_ids[sub_pos] != 0 :
                print('Error -- this sub position is already filled! (line EX)')
            element_ids[sub_pos] = 2000 + lid

        # Now loads
        load_ids = np.where(env.load_to_subid == sub)[0]
        for lid in load_ids :
            sub_pos = env.load_to_sub_pos[lid]
            if element_ids[sub_pos] != 0 :
                print('Error -- this sub position is already filled! (loads)')
            element_ids[sub_pos] = 3000 + lid

        # Now generators
        gen_ids = np.where(env.gen_to_subid == sub)[0]
        for gid in gen_ids :
            sub_pos = env.gen_to_sub_pos[gid]
            if element_ids[sub_pos] != 0 :
                print('Error -- this sub position is already filled! (gens)')
            element_ids[sub_pos] = 4000 + gid

        sub_classes.append(Substation(sub,env.sub_info[sub],element_ids))

    return sub_classes
