load('our-algs.sage')
load('eval.sage')
from ctool import OpCount


#k, A, p, l = 1, 52, 131, 19
k, A, p, l = [3, 8, 31, 19]
#k, A, p, l = [9, 15, 23, 19]
#k, A, p, l = [1, 25, 79, 23]
#k, A, p, l = [11, 57, 101, 23]
#k, A, p, l = [1, 3, 43, 13]
#k, A, p, l = [3, 7, 23, 13]


F = GF(p)
K.<x> = F.extension(k) # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])
N = E.order()
Eg = EllipticCurve(F, [0, A, 0, 1, 0])

G =  point_finding(A,p,l,k, N)
print(G)
G = [G[0], G[2]]
print("Point G: {}".format(G))
Q = [Eg.random_point()]#generate_points(E, l-1)
print(Q)
print("Points to eval: {}".format(Q))
OpCount.clean()
images, ker = algorithm_1(G, Q, A, l)
print("algorithm_1 images: {}".format(images))
OpCount.print_results()
images = algorithm_2(G, Q, A, l)
print("algorithm_2 images: {}".format(images))
images = algorithm_1_using_alg3(G, Q, A, l)
print("algorithm_1_using_alg3 images: {}".format(images))

images, kernel = algorithm_1_using_alg4(G, Q, A, l)

print("algorithm_1_using_alg4 images: {}".format(images))
OpCount.clean()

K_t = GF(l)
r = K_t(2).multiplicative_order()
r_p = r
if r % 2 == 0:
    r_p = r/2
k_p = k
if k % 2 == 0:
    k_p = k/2

m_l = (l-1)/2
c_f = m_l/k_p
a_lk = (r_p)/(gcd(r_p, k_p))
b_lk = c_f /a_lk

images =algorithm_7(G, l, k, Q, A, a_lk, b_lk, p)
print("evaluate_from_G_norm images: {}".format(images))
OpCount.print_results()
