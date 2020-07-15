from math import sin, cos, pi

import numpy as np
from MathematicModels import Quaternion
from scipy.integrate import quad

ez = 0.1
T = 0.1
N = 1
# lam0 = Quaternion([1, 0, 0, 0])
# lam0 = Quaternion([0.8, 0, 0.6, 0])
lam0 = Quaternion([1, 2, 3, 4]) * Quaternion([1/(30.0)**0.5, 0, 0, 0])
Nb = 0.35
    

def r(phi):
    return 1.0 / (1.0 + ez * cos(phi))

def Nk(phi, k):
    # return phi ** k
    return sin((pi * k * phi) / (2 * T))
    
def dNk(phi, k):
    # return k * (phi ** (k-1))
    return ((pi * k) / (2 * T)) * cos((pi * k * phi) / (2 * T))

def lam(a, phi):
    ans = lam0
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
    ans = lam0 * omega * Quaternion([Nk(phi, s+1), 0, 0, 0])
    # print(ans[idx])
    return ans[idx]
    
def Kvect(phi, s, k, idx):
    assert(0<=idx<4)
    # print("s = ", s, "idx = ", idx)
    
    # k -> k+1 the same mistake as in circle0
    # N=1, ||lam(T)||: 0.075 -> 0.997
    ans = Quaternion([2 * dNk(phi, k+1), 0, 0, 0])
    
    r3 = r(phi) ** 3.0
    omega = Quaternion([0, Nb * r3, 0, 1.0])
    
    ans = ans - omega * Quaternion([Nk(phi, k+1), 0, 0, 0])
    ans = ans * Quaternion([Nk(phi, s+1), 0, 0, 0])
    # print(ans[idx])
    return ans[idx]
    
if __name__ == "__main__":
    print("Galerkin >>>")    
    print(r(3.14))
    
    print(lam([Quaternion([1, 0, 0, 0]), Quaternion([0, 1, 1, 1])], pi/2))
    
    K = [[Quaternion([0, 0, 0, 0]) for _ in range(N)] for _ in range(N)]
    #print("K =", K)
    
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
        f[s] = Quaternion([ quad(fvect, 0, T,args=(s, idx))[0] for idx in range(4) ] )
        print("f[", s, "] =", f[s])
        print("f[", s, "][0] =", f[s][0])
        # 
    
    # lam0[0] = 9
    # print(lam0)
    
    assert(N == 1)
    
    a = K[0][0].getInv()
    print("Inv check (1, 0, 0, 0) =", a * K[0][0])
    
    a = f[0] * a 
    print("a =", a)
    
    print("lam0 = ", lam0, lam0.getNorm())
    print("lam(0) = ", lam([a], 0))
    print("lam(T) = ", lam([a], T), lam([a], T).getNorm())
    
    print("Galerkin <<<")    


