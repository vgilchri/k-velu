from ctool import OpCount
load('aux.sage')
load('our-algs.sage')


def evaluate_from_G(p,k,G,A,l,P, k_points):
    '''
    Returns the evaluation at P of the Kernel polynomial generated by G
     ****** WE are assuming k odd ****
     kernel_points is the S0
     '''
    Gx = G[0]
    Px = P[0]
    #k_points = algorithm_4(G, l, m, a, A)
    res = Px - k_points[0][0]
    OpCount.op("add", str(k))
    for i in range (1,len(k_points)) : #Multiplies all generators of Galois orbits
        res=res*(Px-k_points[i][0])
        OpCount.op("mult", str(k))
        OpCount.op("add", str(k))
    power = 1
    for i in range(1, k):
        power = power + pow(p, i)
        OpCount.op("frob", str(k))
    res = pow(res, power)
    res = K(res)
    return res
