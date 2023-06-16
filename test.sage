load('our-algs.sage')
load('eval.sage')
from ctool import OpCount


#k, A, p, l = 1, 52, 131, 19
k, A, p, l = [3, 7, 23, 13]



K.<x> = GF(p^k) # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])

G =  point_finding(A,p,l,k)
G = [G[0], G[2]]
print("Point G: {}".format(G))
Q = [E.point([22*x^2 + 14*x + 21 , 6*x^2 + 11*x + 2 , 1])]#generate_points(E, l-1)
print("Points to eval: {}".format(Q))
images = algorithm_1(G, Q, A, l)
print("algorithm_1 images: {}".format(images))
images = algorithm_2(G, Q, A, l)
print("algorithm_2 images: {}".format(images))
images = algorithm_1_using_alg3(G, Q, A, l)
print("algorithm_1_using_alg3 images: {}".format(images))

images, kernel = algorithm_1_using_alg4(G, Q, A, l)

print("algorithm_1_using_alg4 images: {}".format(images))

OpCount.clean()
m = (l-1)/2
a = 1
algorithm_4(G, l, m, a, A)
OpCount.print_results()

OpCount.clean()
images = evaluate_from_G_norm(p,k,G,A,l, Q, kernel)
print("evaluate_from_G_norm images: {}".format(images))
OpCount.print_results()
