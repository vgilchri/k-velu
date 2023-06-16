load('our-algs.sage')
load('eval.sage')
from ctool import OpCount


k, A, p, l = 1, 52, 131, 19
#k, A, p, l = [9, 15, 23, 19]
#k, A, p, l = [11, 10, 29, 23]
#k, A, p, l = [3, 7, 23, 13]
#k, A, p, l = [11, 57, 101, 23]



K.<x> = GF(p^k) # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])

G =  point_finding(A,p,l,k)
G = [G[0], G[2]]
print("Point G: {}".format(G))
Q = [E.random_point()]#generate_points(E, l-1)
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
t = K(2)

m = (l-1)/2
a = 4
algorithm_4(G, l, m, a, A)
OpCount.print_results()

images = evaluate_from_G_norm(p,k,G,A,l, Q, kernel)
print("evaluate_from_G_norm images: {}".format(images))
OpCount.print_results()
