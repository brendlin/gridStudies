
import CommonHelpers
import BusTopologyHelpers
import numpy as np
import itertools

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

        # Precompute the valid boolean bus states, for all combos of lines on/off and bus configs.
        _lineOffCombos = []
        for i in range(len(self.LocalLineIndices())+1) :
            _lineOffCombos += itertools.combinations(self.LocalLineIndices(),i)

        all_bus_combos = reversed(range(CommonHelpers.FullyConnectedBitset(self.nElements)+1))

        for sub_bitset in all_bus_combos :

            # If it is not a valid bus state at all (regardless of disconnected lines),
            # then do not put the bitset into the validityCache keys at all.
            # We will use the validityCache to get the list of valid possible bus states.
            if not BusTopologyHelpers.IsValidBooleanBusState(sub_bitset,
                                                             self.nElements,
                                                             i_gens = self.LocalGeneratorIndices(),
                                                             i_loads = self.LocalLoadIndices(),
                                                             items_disconnected = []) :
                continue

            self.validityCache[sub_bitset] = dict()

            for _lineOffTuple in _lineOffCombos :
                _lineOffList = list(_lineOffTuple)
                v = BusTopologyHelpers.IsValidBooleanBusState(sub_bitset,
                                                              self.nElements,
                                                              i_gens = self.LocalGeneratorIndices(),
                                                              i_loads = self.LocalLoadIndices(),
                                                              items_disconnected = _lineOffList,
                                                              verbose=False)
                lineOffKey = 0
                for i in _lineOffList :
                    lineOffKey += (0b1 << i)
                self.validityCache[sub_bitset][lineOffKey] = v
                #print('LocalLineIndices:',self.LocalLineIndices(),'checking',_lineOffList,valid)

        return

    def GetValidBusStates(self,lineOnBits=-1) :
        # Get the valid bus states (ignoring effects of turning off lines)
        # These should be filled in at construction.
        if lineOnBits < 0 :
            return list(self.validityCache.keys())

        valid_states = []
        tmp_saveBusConfig = self.currentBusConfig
        for i in self.validityCache.keys() :
            self.SetBusConfig(i)
            if self.IsValidBooleanBusState(lineOnBits=lineOnBits) :
                valid_states.append(i)

        self.currentBusConfig = tmp_saveBusConfig
        return valid_states

    def printValidityCache(self) :
        total = 0
        for kbus in self.validityCache.keys() :
            kbus_str = "0b{:0{}b}".format(kbus,self.nElements)
            invalid = []
            valid = []
            for kline in self.validityCache[kbus].keys() :
                kline_str = "0b{:0{}b}".format(kline,len(self.LocalLineIndices()))
                total += 1
                if self.validityCache[kbus][kline] :
                    valid.append(kline_str)
                else :
                    invalid.append(kline_str)
            print('Bus:',kbus_str)
            print(' - Valid (lineOff) ({})   : {}'.format(len(valid),valid))
            print(' - Invalid: (lineOff) ({}): {}'.format(len(invalid),invalid))
        print('Total states stored in the cache for this bus:',total)
        return

    # "Local" index here refers to the local bus bits here
    def LocalGeneratorIndices(self) :
        return list(np.where(self.elementIDs//1000 == 4)[0])
        
    def LocalLoadIndices(self) :
        return list(np.where(self.elementIDs//1000 == 3)[0])

    def LocalLineIndices(self) :
        return list(np.where(self.elementIDs//1000 <= 2)[0])

    def GetLineID(self,elementID) :
        # For line IDs of the form 2005.001 return "5" (the lineID)
        return int(elementID%1000)
        
    def SetBusConfig(self,bits) :
        self.currentBusConfig = bits
        #if bits not in self.validityCache.keys() :
        #    self.validityCache[bits] = dict()
        return

    def ApplyBusConfig(self,bits,adjacency_matrix_class,verbose=False) :

        tmp_saveBusConfig = self.currentBusConfig
        self.SetBusConfig(bits)
        if not self.IsValidBooleanBusState(lineOnBits=adjacency_matrix_class.lineOnBits) :
            tmp = 'Warning: tried to apply an invalid bus state 0b{:0{}b}. Doing nothing.'
            print(tmp.format(bits,self.nElements))
            self.currentBusConfig = tmp_saveBusConfig
            return

        if verbose :
            print('Switching from 0b{:0{}b} to 0b{:0{}b}'.format(self.currentBusConfig,
                                                                 self.nElements,
                                                                 bits,
                                                                 self.nElements))
        toBus1 = []
        toBus2 = []
        for i in range(self.nElements) :
            if (bits & (0b1 << i)) == (0b1 << i) :
                toBus1.append(self.elementIDs[i])
            else :
                toBus2.append(self.elementIDs[i])

        if verbose :
            print('On bus 1:',toBus1)
            print('On bus 2:',toBus2)

        self.currentBusConfig = bits

        if len(toBus1) :
            adjacency_matrix_class.SetListOfElementsToBusN(self.index,1,toBus1)
        if len(toBus2) :
            adjacency_matrix_class.SetListOfElementsToBusN(self.index,2,toBus2)

        return

    def IsValidBooleanBusState(self,lineOnBits=-1) :

        # If it is not in the validity cache, then it is not a valid bus state, period.
        #print('Checking if {} is in'.format(self.currentBusConfig),self.validityCache.keys())
        if self.currentBusConfig not in self.validityCache.keys() :
            #print('It is not! Getting out of here!')
            return False

        localDisconnIndices = []
        cacheKey = 0

        if lineOnBits > 0 :
            for li in self.LocalLineIndices() :
                lineID = self.GetLineID(self.elementIDs[li])
                isOn = lineOnBits & (0b1 << lineID)
                #print('{} is {}'.format(lineID,'on' if isOn else 'off'))
                if isOn :
                    continue
                else :
                    cacheKey += (0b1 << li)
                    localDisconnIndices.append(li)

        # Otherwise, all results are now precomputed.
        # print('Looking for bus 0b{:0{}b}, key 0b{:0{}b}'.format(self.currentBusConfig,
        #                                                         self.nElements,
        #                                                         cacheKey,
        #                                                         len(self.LocalLineIndices())))
        return self.validityCache[self.currentBusConfig][cacheKey]


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
            links_to_sub = env.line_ex_to_subid[lid]
            element_ids[sub_pos] = 1000 + lid + np.round(links_to_sub / 1000.,3)

        # Find the line (EX) ids that link to this sub
        line_ex_ids = np.where(env.line_ex_to_subid == sub)[0]
        for lid in line_ex_ids :
            sub_pos = env.line_ex_to_sub_pos[lid]
            if element_ids[sub_pos] != 0 :
                print('Error -- this sub position is already filled! (line EX)')
            links_to_sub = env.line_or_to_subid[lid]
            element_ids[sub_pos] = 2000 + lid + np.round(links_to_sub / 1000.,3)

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
