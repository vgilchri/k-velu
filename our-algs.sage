load('aux.sage')

def normalize_images(images, K):
    normlized = []
    for image in images:
        r1 = K(image[0])/K(image[1])
        r2 = K(1)
        normlized.append([r1, r2])
    return normlized


def kernel_points(P, A, d):
	# alg 2, p. 11
	# given a generator of the kernel, P, and the montgomery curve constant of the domain curve,
	# returns the first d multiples of P
    kernel = [P]
    if d >= 2:
        kernel.append(xDBL(P,A))
        for i in range(2,d):
            temp = xADD(kernel[i-1],P,kernel[i-2])
            kernel.append(temp)
    # if d==2, then kernel will have two elts: P which is added at the start, and infinity which is added at the end
    return kernel

def algorithm_3(G, A, l):
    kernel = [G]
    double_G = xDBL(G, A)
    kernel.append(double_G)
    for i in range(2, (l-1)/2):
        n_G = xADD(kernel[i-1], G, kernel[i-2])
        kernel.append(n_G)
    return kernel

#m = (l-1)/2
#a = #A_l where A_l = <2>
def algorithm_4(G, l, m, a, A):
    ''' We assume that l satisfies the condition'''
    b = int(m//a)
    kernel = []
    for i in range(b):
        idx = (a*i)
        if i == 0:
            kernel.append(G)
        else:
            T = kernel[idx]
            T1 = kernel[idx-1]
            P = xADD(T, T1, T1)
            kernel.append(P)
        for j in range(1, a):
            tt = (idx+j)-1
            T = kernel[tt]
            P = xDBL(T, A)
            kernel.append(P)
    return kernel


def algorithm_1(G, eval_points, A, l):
    hat_points = []
    K = G[0].parent()
    kernel = kernel_points(G, A, (l-1)/2)
    #print(kernel)
    for i in range(0, (l-1)/2):
        x_hat = kernel[i][0] + kernel[i][1]
        z_hat = kernel[i][0] - kernel[i][1]
        OpCount.op("A", str(k))
        OpCount.op("A", str(k))
        hat_points.append([x_hat, z_hat])

    images = []
    for i in range(len(eval_points)):
        u_hat = eval_points[i][0] + eval_points[i][1]
        v_hat = eval_points[i][0] - eval_points[i][1]
        OpCount.op("A", str(k))
        OpCount.op("A", str(k))
        u_prime = K(1)
        v_prime = K(1)
        for j in range(0, (l-1)/2):
            t0, t1 = criss_cross(hat_points[j][0], hat_points[j][1], u_hat, v_hat)
            u_prime = t0*u_prime
            v_prime = t1*v_prime
            OpCount.op("M", str(k))
            OpCount.op("M", str(k))
        u_prime = eval_points[i][0] * (u_prime**2)
        v_prime = eval_points[i][1] * (v_prime**2)
        OpCount.op("C", str(k))
        OpCount.op("C", str(k))
        OpCount.op("S", str(k))
        OpCount.op("S", str(k))
        images.append([u_prime, v_prime])

    images = normalize_images(images, K)

    return images, kernel

def algorithm_2(G, eval_points, A, l):
    hat_points = []
    K = G[0].parent()
    kernel = kernel_points(G, A, (l-1)/2)
    n = len(eval_points)

    tmp_eval = []
    tmp_eval_prime = []
    for i in range(n):
        U_hat = eval_points[i][0]+eval_points[i][1]
        V_hat = eval_points[i][0]-eval_points[i][1]
        tmp_eval.append([U_hat, V_hat])
        tmp_eval_prime.append([K(1), K(1)])

    for P in kernel:
        x_hat = P[0] + P[1]
        z_hat = P[0] - P[1]
        for i in range(n):
            t0, t1 = criss_cross(x_hat, z_hat, tmp_eval[i][0], tmp_eval[i][1])
            u_i = tmp_eval_prime[i][0]*t0
            v_i = tmp_eval_prime[i][1]*t1
            tmp_eval_prime[i] = [u_i, v_i]

    images = []
    for i in range(len(eval_points)):
        ui_prime = tmp_eval_prime[i][0]**2
        vi_prime = tmp_eval_prime[i][1]**2
        u_prime = eval_points[i][0] * ui_prime
        v_prime = eval_points[i][1] * vi_prime
        images.append([u_prime, v_prime])
    print(images)
    images = normalize_images(images, K)

    return images

def algorithm_1_using_alg3(G, eval_points, A, l):
    hat_points = []
    K = G[0].parent()
    kernel = algorithm_3(G, A, l)
    #print("Kernel alg 3:{} ".format(kernel))
    for i in range(0, (l-1)/2):
        x_hat = kernel[i][0] + kernel[i][1]
        z_hat = kernel[i][0] - kernel[i][1]
        hat_points.append([x_hat, z_hat])

    images = []
    for i in range(len(eval_points)):
        u_hat = eval_points[i][0] + eval_points[i][1]
        v_hat = eval_points[i][0] - eval_points[i][1]
        OpCount.op("add", str(k))
        OpCount.op("add", str(k))
        u_prime = K(1)
        v_prime = K(1)
        for j in range(0, (l-1)/2):
            t0, t1 = criss_cross(hat_points[j][0], hat_points[j][1], u_hat, v_hat)
            u_prime = t0*u_prime
            v_prime = t1*v_prime
            OpCount.op("mult", str(k))
            OpCount.op("mult", str(k))
        u_prime = eval_points[i][0] * (u_prime**2)
        v_prime = eval_points[i][1] * (v_prime**2)
        OpCount.op("mult", str(k))
        OpCount.op("mult", str(k))
        OpCount.op("square", str(k))
        OpCount.op("square", str(k))
        images.append([u_prime, v_prime])
    images = normalize_images(images, K)

    return images

def algorithm_1_using_alg4(G, eval_points, A, l):
    hat_points = []
    K = G[0].parent()
    m = (l-1)/2
    a = m
    kernel = algorithm_4(G, l, m, a, A)
    #print("Kernel alg 4:{} ".format(kernel))
    for i in range(0, (l-1)/2):
        x_hat = kernel[i][0] + kernel[i][1]
        z_hat = kernel[i][0] - kernel[i][1]
        OpCount.op("A", str(k))
        OpCount.op("A", str(k))
        hat_points.append([x_hat, z_hat])

    images = []
    for i in range(len(eval_points)):
        u_hat = eval_points[i][0] + eval_points[i][1]
        v_hat = eval_points[i][0] - eval_points[i][1]
        OpCount.op("A", str(k))
        OpCount.op("A", str(k))
        u_prime = K(1)
        v_prime = K(1)
        for j in range(0, (l-1)/2):
            t0, t1 = criss_cross(hat_points[j][0], hat_points[j][1], u_hat, v_hat)
            u_prime = t0*u_prime
            v_prime = t1*v_prime
            OpCount.op("M", str(k))
            OpCount.op("M", str(k))
        u_prime = eval_points[i][0] * (u_prime**2)
        v_prime = eval_points[i][1] * (v_prime**2)
        OpCount.op("C", str(k))
        OpCount.op("C", str(k))
        OpCount.op("S", str(k))
        OOpCount.op("S", str(k))
        images.append([u_prime, v_prime])
    images = normalize_images(images, K)

    return images, kernel

def algorithm_6(G, eval_points, A, l):
    hat_points = []
    K = G[0].parent()
    m = (l-1)/2
    a = m
    kernel = algorithm_4(G, l, m, a, A)
    #print("Kernel alg 4:{} ".format(kernel))
    for i in range(0, (l-1)/2):
        x_hat = kernel[i][0] + kernel[i][1]
        z_hat = kernel[i][0] - kernel[i][1]
        OpCount.op("add", str(k))
        OpCount.op("add", str(k))
        hat_points.append([x_hat, z_hat])

    images = []
    for i in range(len(eval_points)):
        u_hat = eval_points[i][0] + eval_points[i][1]
        v_hat = eval_points[i][0] - eval_points[i][1]
        OpCount.op("add", str(k))
        OpCount.op("add", str(k))
        u_prime = K(1)
        v_prime = K(1)
        for j in range(0, (l-1)/2):
            t0, t1 = criss_cross(hat_points[j][0], hat_points[j][1], u_hat, v_hat)
            u_prime = t0*u_prime
            v_prime = t1*v_prime
            OpCount.op("mult", str(k))
            OpCount.op("mult", str(k))
        u_prime = eval_points[i][0] * (u_prime**2)
        v_prime = eval_points[i][1] * (v_prime**2)
        OpCount.op("mult", str(k))
        OpCount.op("mult", str(k))
        OpCount.op("square", str(k))
        OpCount.op("square", str(k))
        images.append([u_prime, v_prime])
    images = normalize_images(images, K)

    return images, kernel

def algorithm_5(G, a_l, b_k, A):
    kern = []
    for i in range(b_k):
        if i == 0:
            kern.append(G)
        else:
            P = xADD(kern[a_l*i], kern[(a_l*i)-1], kern[(a_l*i)-1])
            kern.append(P)
        for j in range(1, a_l):

            T = xDBL(kern[((a_l*i)+j)-1], A)
            kern.append(T)
    return kern

def c_norm(alpha, k, q):
    t = copy(alpha)
    for i in range(k-1):
        t = (t**q) * alpha
        OpCount.op("M", str(k))
        OpCount.op("F", str(k))

    return t

def algorithm_7(G, l, k, eval_points, A, a_l, b_k, p):
    K = G[0].parent()

    kern_point = algorithm_5(G, a_l, b_k, A)

    hat_points = []
    for P in kern_point:
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

    normalize_images(images, K)
    return images
