from ctool import OpCount
def xADD(P,Q,R):
	# montgomery xADD
	#OpCount.op("xADD", str(k))
	xP,zP = P
	xQ,zQ = Q
	xR,zR = R
	U = (xP-zP)*(xQ+zQ)
	V = (xP+zP)*(xQ-zQ)
	res1 = zR*((U+V)**2)
	res2 = xR*((U-V)**2)
	OpCount.op("mult", str(k))
	OpCount.op("mult", str(k))
	OpCount.op("mult", str(k))
	OpCount.op("mult", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("square", str(k))
	OpCount.op("square", str(k))
	if res2 == 0:
		res1 = 1
		res2 = 0
	else:
		res1 = res1/res2
		res2 = 1
	return [res1, res2]


def xDBL(P,A):
	# montgomery xDBL
	xP,zP = P
	R = (xP+zP)**2
	S = (xP-zP)**2
	t = xP*zP # T/4 from s.s. paper
	r1 = R*S
	r2 = 4*t*(S+(A+2)*t)
	OpCount.op("mult", str(k))
	OpCount.op("mult", str(k))
	OpCount.op("square", str(k))
	OpCount.op("square", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("C", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	if r2 == 0:
		r1 = 1
		r2 = 0
	else:
		r1 = r1/r2
		r2 = 1
	return [r1, r2]

def criss_cross(a,b,c,d):
	# alg 1, p. 11
	# performs a small computation on the inputs
	# cost: 2M + 2a
	t1 = a*d
	OpCount.op("mult", str(k))
	t2 = b*c
	OpCount.op("mult", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	return (t1+t2, t1-t2)

#def kernel_points(P, A, d):
	# alg 2, p. 11
	# given a generator of the kernel, P, and the montgomery curve constant of the domain curve,
	# returns the first d multiples of P
#	kernel = [P]
#	if d >= 2:
		#kernel.append(xDBL(P,A))
		#for i in range(2,d):
		#	temp = xADD(kernel[i-1],P,kernel[i-2])
		#	kernel.append(temp)
	# if d==2, then kernel will have two elts: P which is added at the start, and infinity which is added at the end
	#return kernel

from ctool import OpCount

# Normal point-finding algo

def point_finding(A,p,l,k):
# checks
    if (l-1 % k == 0):
        return "k ne divise pas l-1"
    if (k > 13):
        return "k est trop grand"
    K = GF(p^k, 'x')
    print("Field: {}".format(K))
    # sample a random point
    E = EllipticCurve(K, [0, A, 0, 1, 0])
    P_l =E.random_point()
    N = E.order()
    if N % l**2 == 0:
        t = N // l**2
        OpCount.op("div", str(k))
    else:
        t = N//l
        OpCount.op("div", str(k))

    while P_l.order() != l:
        P = E.random_point()
        OpCount.op("mult", str(k))
        P_l = t*P
    return P_l

def generate_points(E, n):
    G = E.gens()[0]
    points =  []
    for i in range(1, n+1):
        P = i*G
        points.append([P[0], P[2]])
    return points

def frob_power(K,E,P,power):
    frob = K.frobenius_endomorphism(power)
    R = E(frob(P[0]), frob(P[1]))
    return R

# returns the image of a point P in E by the endo delta

# delta is defined as the division of (X^k - 1) and the kth cyclotomic polynomial

def delta(P,E,k,p):
    R = "error, wrong k"
    K = GF(p**k, 'x')
    if k == 1:
        return P
    if (k in {2,3,5,7,11}):
        R = frob_power(K,E,P,1) - P
    elif k==4 :
        R = frob_power(K,E,P,2) - P
    elif k==6 :
        S = frob_power(K,E,P,3) - P
        R = frob_power(K,E,S,1) + S
    elif k==8 :
        R = frob_power(K,E,P,4) - P
    elif k==9 :
        R = frob_power(K,E,P,3) - P
    elif k==10 :
        S = frob_power(K,E,P,5) - P
        R = frob_power(K,E,S,1) + S
    elif k==12 :
        S = frob_power(K,E,P,6) - P
        R = frob_power(K,E,S,2) + S
    return R

def get_h_k(A,p,k) :

    K = GF(p**k)

    E = EllipticCurve(K, [0, A, 0, 1, 0])

    Nk = E.order()
    N = Nk

    if (k in {2,3,5,7,11}) :

        K1 = GF(p)

        E1 = EllipticCurve(K1, [0, A, 0, 1, 0])

        N1 = E1.order()

        N= Nk//N1


    elif k==4 :

        K2 = GF(p**2)

        E2 = EllipticCurve(K2, [0, A, 0, 1, 0])

        N2 = E2.order()

        N= Nk//N2


    elif k==6 :

        K1 = GF(p)

        E1 = EllipticCurve(K1, [0, A, 0, 1, 0])

        N1 = E1.order()

        K2 = GF(p**2)

        E2 = EllipticCurve(K2, [0, A, 0, 1, 0])

        N2 = E2.order()

        K3 = GF(p**3)

        E3 = EllipticCurve(K3, [0, A, 0, 1, 0])

        N3 = E3.order()

        N= (Nk*N1)//(N2*N3)


    elif k==8 :

        K4 = GF(p**4)

        E4 = EllipticCurve(K4, [0, A, 0, 1, 0])

        N4 = E4.order()

        N= Nk//N4


    elif k==9 :

        K3 = GF(p**3)

        E3 = EllipticCurve(K3, [0, A, 0, 1, 0])

        N3 = E3.order()

        N= Nk//N3


    elif k==10 :

        K1 = GF(p)

        E1 = EllipticCurve(K1, [0, A, 0, 1, 0])

        N1 = E1.order()

        K2 = GF(p**2)

        E2 = EllipticCurve(K2, [0, A, 0, 1, 0])

        N2 = E2.order()

        K5 = GF(p**5)

        E5 = EllipticCurve(K5, [0, A, 0, 1, 0])

        N5 = E5.order()

        N= (Nk*N1)//(N2*N5)

    elif k==12 :

        K2 = GF(p**2)

        E2 = EllipticCurve(K2, [0, A, 0, 1, 0])

        N2 = E2.order()

        K4 = GF(p**4)

        E4 = EllipticCurve(K4, [0, A, 0, 1, 0])

        N4 = E4.order()

        K6 = GF(p**6)

        E6 = EllipticCurve(K6, [0, A, 0, 1, 0])

        N6 = E6.order()

        N= (Nk*N2)//(N4*N6)


    return N

def optimized_point_finding(A,p,k,l, K):
	E = EllipticCurve(K, [0, A, 0, 1, 0])
	P_l = E(0)
	N = get_h_k(A,p,k)
	print("N: {}".format(N))
	if N % l**2 == 0:
		t = N // l**2
	else:
		t = N//l
	while P_l.order() != l:
		P = E.random_point()
		P = delta(P,E,k,p)
		P_l = t*P
		print("P: {}".format(P_l))

	return P_l
