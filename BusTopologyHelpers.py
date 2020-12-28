
# Module for imposing bus rules, given a bitset representation.

import CommonHelpers

def RemoveBitFromFlag(flag,bit_to_remove,verbose=False) :
    # Remove a bit (count starting from 0)
    # Returns a flag of length n-1

    # If the bit value is higher than the flag length, then do nothing.
    if (0b1 << bit_to_remove) > flag :
        return flag
    
    bits_above = flag >> (bit_to_remove+1) << (bit_to_remove+1)
    bits_below = flag - bits_above - (0b1 << bit_to_remove)
    new_bits = (bits_above >> 1) + bits_below
    if verbose :
        print('The bits  : 0b{:06b}'.format(flag))
        print('Removing  : 0b{:06b}'.format(0b1 << bit_to_remove))
        print('Bits below: 0b{:06b}'.format(bits_below))
        print('Bits above: 0b{:06b}'.format(bits_above))
        print('New bits  : 0b{:06b}'.format(new_bits))
    return new_bits


def IsValidBooleanBusState(flag,nbits,i_gens=[],i_loads=[],items_disconnected=[],verbose=False) :
    # Bus rules are imposed to make sure they are:
    # (a) unique (i.e. not duplicate)
    # (b) legal (i.e. there must be at least 2 items on a bus)
    # (c) unique and consistent in the event of a disconnected line
    #  - To impose (c), we require that a disconnected line be on bus "1" instead of bus "0".
    #    In the case a disconnected line is on bus "0", we consider it a duplicate option.
    #    This is just a simplifying convention.
    
    flag_str = "{:0{}b}".format(flag,nbits)
    #print(flag_str)
    n_ones = flag_str.count('1')
    n_zeros = flag_str.count('0')

    # By convention, the leading item must be on bus "1" (to remove 1 <--> 0 symmetry duplicates)
    if not (flag & (0b1 << (nbits-1))) :
        if verbose : print('Excluding {:0{}b} by symmetry convention - Flag does not have a leading 1.'.format(flag,nbits))
        return False

    # You cannot have a single line on a bus
    if ((n_ones) == 1 or (n_zeros == 1)) and nbits > 1:
        if verbose : print('Excluding {:0{}b} by illegal bus manoeuvre - {:d} item(s) on bus 1 and {:d} item(s) on bus 0.'.format(flag,nbits,n_ones,n_zeros))
        return False

    # Test to make sure generators or loads (externals) are not islanded.
    # (In other words, does one bus only contain externals and disconnected lines?)
    if len(i_gens + i_loads + items_disconnected) > 0 :

        externals_and_disconnecteds = 0
        for i in i_gens + i_loads + items_disconnected :
            externals_and_disconnecteds += 0b1 << i

        if flag and (flag & externals_and_disconnecteds == flag) :
            if verbose : print('Excluding {:0{}b} because of islanded externals on bus 1.'.format(flag,nbits))
            return False

        flag_inverted = CommonHelpers.FullyConnectedBitset(nbits) - flag
        if flag_inverted and (flag_inverted & externals_and_disconnecteds == flag_inverted) :
            if verbose : print('Here: Excluding {:0{}b} because of islanded externals on bus 0.'.format(flag,nbits))
            return False

    # If lines are disconnected, then by convention they must be set to bus "1" here
    reduced_flag = flag
    for i,item in enumerate(sorted(items_disconnected,reverse=True)) :

        if not (flag & (0b1 << item)) :
            if verbose : print('Excluding {:0{}b} by symmetry convention - disconnected line {} is not on bus 1'.format(flag,nbits,item))
            return False
        
        # Delete this bit and check that the remaining bits are valid
        if verbose : print('Testing old flag {:0b}'.format(flag))
        reduced_flag = RemoveBitFromFlag(reduced_flag,item,verbose=verbose)
        if verbose : print('Old flag was: {:0b} new flag is: {:0b}'.format(flag,reduced_flag))
        if (not IsValidBooleanBusState(reduced_flag,nbits-1-i)) :
            return False
    
    return True
