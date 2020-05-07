from math import sin, cos, pi
import numpy as np
from MathematicModels import Quaternion
from Gauss import solveLinear
from scipy.integrate import quad, odeint, ode
import matplotlib.pyplot as plt

ez = 0.01
T = 4*pi
N = 1
# lam0 = Quaternion([1, 0, 0, 0])
# lam0 = Quaternion([0.8, 0, 0.6, 0])
lam0 = Quaternion([1, 2, 3, 4]) * Quaternion([1/(30.0)**0.5, 0, 0, 0])
Nb = 0.35
    
def lamCircle(phi):
    omega = Quaternion([0, Nb, 0, 1.0])
    om = (1 + Nb * Nb) ** 0.5
    
    omega = omega * Quaternion([sin(0.5 * om * phi)/om, 0, 0, 0])
    ans = Quaternion([cos(0.5 * om * phi), 0, 0, 0])
    ans = ans + omega
    ans = lam0 * ans
    return ans
    
def dlamCircle(phi):
    omega = Quaternion([0, Nb, 0, 1.0])
    om = (1 + Nb * Nb) ** 0.5
    omega = omega * Quaternion([0.5*cos(0.5 * om * phi), 0, 0, 0])
    ans = Quaternion([-0.5 * om * sin(0.5 * om * phi), 0, 0, 0])
    ans = ans + omega
    ans = lam0 * ans
    return ans
    

def r(phi):
    return 1.0 / (1.0 + ez * cos(phi))

def Nk(phi, k):
    # return phi ** k
    return sin((pi * k * phi) / (2 * T))
    # return sin((k * phi) / (2 * T))
    
def dNk(phi, k):
    # return k * (phi ** (k-1))
    return ((pi * k) / (2 * T)) * cos((pi * k * phi) / (2 * T))
    # return ((k) / (2 * T)) * cos((k * phi) / (2 * T))

def lam(a, phi):
    ans = lamCircle(phi)
    for k in range(len(a)):
        # print("ak", a[k])
        # Nkq = Quaternion([Nk(phi, k+1), 0, 0, 0])
        # print("Nkq", Nkq)
        # mul = a[k] * Nkq
        # print("mul", mul)
        ans = ans + a[k] * Quaternion([Nk(phi, k+1), 0, 0, 0])
    return ans

def fvect(phi, s, idx):
    assert(0<=idx<4)
    # print("s = ", s, "idx = ", idx)
    r3 = r(phi) ** 3.0
    omega = Quaternion([0, Nb * r3, 0, 1.0])
    ans = lamCircle(phi) * omega - Quaternion([2.0, 0, 0, 0]) * dlamCircle(phi)
    
    ans = ans * Quaternion([Nk(phi, s+1), 0, 0, 0])
    # print(ans[idx])
    return ans[idx]
    
def Kvect(phi, s, k, idx):
    assert(0<=idx<4)
    # print("s = ", s, "idx = ", idx)
    
    ans = Quaternion([2 * dNk(phi, k), 0, 0, 0])
    
    r3 = r(phi) ** 3.0
    omega = Quaternion([0, Nb * r3, 0, 1.0])
    
    ans = ans - omega * Quaternion([Nk(phi, k+1), 0, 0, 0])
    ans = ans * Quaternion([Nk(phi, s+1), 0, 0, 0])
    # print(ans[idx])
    return ans[idx]
    
# Cauchy
def eqEllipse(lam, phi):
    #print("lam =", lam)
    lam = Quaternion(lam)
    r3 = r(phi) ** 3.0
    omega = Quaternion([0, Nb * r3, 0, 1.0])
    ans = lam * omega
    ans = ans * Quaternion([0.5, 0, 0, 0])
    return [ans[idx] for idx in range(4)]
        
    
if __name__ == "__main__":
    print("Galerkin >>>")    
    print(r(3.14))
    
    print(lam([Quaternion([1, 0, 0, 0]), Quaternion([0, 1, 1, 1])], pi/2))
    
    K = [[Quaternion([0, 0, 0, 0]) for _ in range(N)] for _ in range(N)]
    print("K =", K)
    
    for s in range(N):
        for k in range(N):
            K[s][k] = Quaternion([ quad(Kvect, 0, T,args=(s, k, idx))[0] for idx in range(4) ] )
            print("K[{0}, {1}] = {2}".format(s, k, K[s][k]))
    
    f = [Quaternion([0, 0, 0, 0]) for _ in range(N)]
    # print("f =", f)
    print("---")
    for s in range(N):
        # res = [ quad(fvect, 0, T,args=(s, idx))[0] for idx in range(4) ]
        # print("res = ", res)
        f[s] = [Quaternion([ quad(fvect, 0, T,args=(s, idx))[0] for idx in range(4) ] )]
        print("f[", s, "] =", f[s])
        print("f[", s, "][0] =", f[s][0])
    
    print("f =", f) 
    
    # lam0[0] = 9
    # print(lam0)
    
    assert(N == 1)
    
    a = K[0][0].getInv()
    print("Inv check (1, 0, 0, 0) =", a * K[0][0])
    
    a = f[0][0] * a 
    # a = Quaternion([0, 0, 0, 0])
    print("a =", a)
    
    aGauss = solveLinear(K, f)
    for idx, elem in enumerate(aGauss):
        print("aGauss[{0}] = {1}".format(idx, list(map(str, elem))))
    
    print("lam0 = ", lam0, lam0.getNorm())
    print("lam(0) = ", lam([a], 0))
    print("lam(T) = ", lam([a], T), lam([a], T).getNorm())
    
    #Cauchy
    # circle [0.031046560091828313, 0.031103992683832533, 0.031161212781716421, 0.031218202529294029, 0.031274944260239294]
    # poly [0.029317544396431784, 0.029377607905505098, 0.029437517515314045, 0.029497254701257254, 0.029556801152425299]
    # sin [0.029878477588461995, 0.02994193916572772, 0.030005212904268913, 0.030068280403747726, 0.030131123479297209]
    print("Cauchy >>>")
    phi = np.linspace(0, T, 1001)
    sol = odeint(eqEllipse, [lam0[idx] for idx in range(4)], phi, rtol=1e-15)
    print(sol[:5], "\n\n", sol[-5:])
    
    lastQ = Quaternion(sol[-1])
    print(lastQ, lastQ.getNorm())
    
    print(lam([a], T), lam([a], T).getNorm())
    
    err = []
    for cur in zip(phi, sol):
        l = lam([a], cur[0])
        diff = Quaternion(cur[1]) - l
        err.append(diff.getNorm() ** .5)
        
    print(err[:5])
    print(err[-5:])
    
    plt.plot(phi, err, label ='err')
    plt.legend(loc='best')
    plt.xlabel('phi')
    plt.grid()
    plt.show()
    
    print("Cauchy <<<")
    
    print("Galerkin <<<")    


