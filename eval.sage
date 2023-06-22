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


def nor(U, V, S0, k):
    Xs = []
    Zs = []
    for i in range(k):
        OpCount.op("frob", str(k))
        OpCount.op("frob", str(k))
        for P in S0:
            Px = pow(P[0], p**i)
            Xs.append(Px)
            Pz = pow(P[1], p**i)
            Zs.append(Pz)
    U_prime = K(1)
    V_prime = K(1)
    for i in range(0, len(Xs)):
        U_prime = U_prime*Xs[i]
        V_prime = V_prime*Zs[i]
    return U_prime, V_prime



def norm_c(alpha, p, k):
    t = copy(alpha)
    for i in range(k-1):
        t = (t**q) * alpha
        OpCount.op("M", str(k))
        OpCount.op("F", str(k))

    return t

def evaluate_from_G_norm(p,k,G,A,l, eval_points, k_points):
    '''
    Returns the evaluation at eval_points of the Kernel points generated by G
     ****** WE are assuming k odd ****
     k_points is the S0
     '''
    hat_points = []
    for P in k_points:
        X_hat = P[0] + P[1]
        Z_hat = P[0] - P[1]
        OpCount.op("A", str(k))
        OpCount.op("A", str(k))
        hat_points.append([X_hat, Z_hat])
    images = []
    for Q in eval_points:
        U_hat = Q[0] + Q[1]
        V_hat = Q[0] - Q[1]
        OpCount.op("A", str(k))
        OpCount.op("A", str(k))
        U_prime = K(1)
        V_prime = K(1)
        for H in hat_points:
            t0,t1 = criss_cross(H[0], H[1], U_hat, V_hat)
            U_prime = t0*U_prime
            V_prime = t1*V_prime
            OpCount.op("M", str(k))
            OpCount.op("M", str(k))
        U_prime = c_norm(U_prime, k, p)
        V_prime = c_norm(V_prime, k, p)
        U_prime = U_prime**2
        V_prime = V_prime**2
        OpCount.op("S", str(k))
        OpCount.op("S", str(k))
        U_prime = Q[0]*U_prime
        OpCount.op("C", str(k))
        V_prime = Q[1]*V_prime
        OpCount.op("C", str(k))
        images.append([U_prime, V_prime])
    images = normalize_images(images, K)
    return images
