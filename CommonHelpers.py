
def FullyConnectedBitset(n_bits) :
    x = 0
    for i in range(n_bits) :
        x += 0b1 << i
    return x
