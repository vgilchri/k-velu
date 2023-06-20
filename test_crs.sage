load('our-algs.sage')
load('eval.sage')
from ctool import OpCount




#k, A, p, l = 1, 52, 131, 19
#k, A, p, l = [3, 8, 31, 19]
#k, A, p, l = [9, 15, 23, 19]
#k, A, p, l = [11, 10, 29, 23]
#k, A, p, l = [11, 57, 101, 23]
#k, A, p, l = [1, 3, 43, 13]
#k, A, p, l = [3, 7, 23, 13]
k = 1
N = 152
p = 12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769
A = 10861338504649280383859950140772947007703646408372831934324660566888732797778932142488253565145603672591944602210571423767689240032829444439469242521864171


# ---------------------	RUN OUR ALGO ----------------------------------------

K.<x> = GF(p^k) # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])


ll = [ 5, 7, 11, 13, 17, 103]
print("Starting.... ")

for l in ll:
    G =  point_finding(A,p,l,k)
    #G = optimized_point_finding(A,p,k,l,K)
    G = [G[0], G[2]]
    print("Point G: {}".format(G))
    Q = [E.random_point()]#generate_points(E, l-1)
    print("Points to eval: {}".format(Q))
    print("Computing l :{}".format(l))
    OpCount.clean()
    images, ker = algorithm_1(G, Q, A, l)
    print("algorithm_1 images: {}".format(images))
    OpCount.print_results()
    #images = algorithm_2(G, Q, A, l)
    #print("algorithm_2 images: {}".format(images))
    #images = algorithm_1_using_alg3(G, Q, A, l)
    #print("algorithm_1_using_alg3 images: {}".format(images))

    #images, kernel = algorithm_1_using_alg4(G, Q, A, l)

    #print("algorithm_1_using_alg4 images: {}".format(images))
    OpCount.clean()
    m = (l-1)/2
    a = m
    kern = algorithm_4(G, l, m, a, A)
    OpCount.print_results()

    images = evaluate_from_G_norm(p,k,G,A,l, Q, kern)
    print("evaluate_from_G_norm images: {}".format(images))
    OpCount.print_results()
