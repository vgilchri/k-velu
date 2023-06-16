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



def get_Norm(k, res):
    power = 1
    for i in range(1, k):
        power = power + pow(p, i)
        OpCount.op("frob", str(k))
    res = pow(res, power)
    res = K(res)
    return res


def evaluate_from_G_norm(p,k,G,A,l, eval_points, k_points):
    '''
    Returns the evaluation at eval_points of the Kernel points generated by G
     ****** WE are assuming k odd ****
     k_points is the S0
     '''
    tmp_eval = []
    tmp_eval_prime = []
    n = len(eval_points)
    for i in range(n):
        U_hat = eval_points[i][0]+eval_points[i][1]
        V_hat = eval_points[i][0]-eval_points[i][1]
        OpCount.op("add", str(k))
        OpCount.op("add", str(k))
        tmp_eval.append([U_hat, V_hat])
        tmp_eval_prime.append([K(1), K(1)])

    for P in k_points:
        x_hat = P[0] + P[1]
        z_hat = P[0] + P[1]
        OpCount.op("add", str(k))
        OpCount.op("add", str(k))
        for i in range(n):
            t0, t1 = criss_cross(x_hat, z_hat, tmp_eval[i][0], tmp_eval[i][1])
            u_i = tmp_eval_prime[i][0]*t0
            OpCount.op("mult", str(k))
            v_i = tmp_eval_prime[i][1]*t1
            OpCount.op("mult", str(k))
            tmp_eval_prime[i] = [u_i, v_i]

    images = []
    for i in range(len(eval_points)):
        ui_prime, vi_prime = nor(eval_points[i][0], eval_points[i][1], k_points, k)
        ui_prime = ui_prime**2
        vi_prime = vi_prime**2
        OpCount.op("square", str(k))
        OpCount.op("square", str(k))
        u_prime = eval_points[i][0] * ui_prime
        v_prime = eval_points[i][1] * vi_prime
        OpCount.op("mult", str(k))
        OpCount.op("mult", str(k))
        images.append([u_prime, v_prime])
    print(images)
    images = normalize_images(images, K)
    return images
