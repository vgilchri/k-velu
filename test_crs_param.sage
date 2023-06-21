load('our-algs.sage')
load('eval.sage')
from ctool import OpCount

k = 3
N = 152
p = 12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769
A = 10861338504649280383859950140772947007703646408372831934324660566888732797778932142488253565145603672591944602210571423767689240032829444439469242521864171

K.<x> = GF(p^k) # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])
l = 19

#G =  point_finding(A,p,l,k)
G = optimized_point_finding(A,p,k,l,K)
G = [G[0], G[2]]
print("Point G: {}".format(G))
Q = [E.random_point()]#generate_points(E, l-1)
print("Points to eval: {}".format(Q))
print("Computing l :{}".format(l))
OpCount.clean()
images, ker = algorithm_1(G, Q, A, l)
print("algorithm_1 images: {}".format(images))
OpCount.print_results()

OpCount.clean()
a_l = 3
b_k = 1
images = algorithm_7(G, l, k, Q, A, a_l, b_k, p)
print("algorithm_7 images: {}".format(images))
OpCount.print_results()
