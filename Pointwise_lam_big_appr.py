"""
    Lambda big eq. No circle.
"""
from math import sin, cos, pi
import numpy as np
from MathematicModels import Quaternion
from Gauss import solveLinear
from copy import deepcopy
from scipy.integrate import quad, odeint, ode
import matplotlib.pyplot as plt

EZ = 0.1
T = pi/2
M = 8
# lam0 = Quaternion([1, 0, 0, 0])
# lam0 = Quaternion([0.8, 0, 0.6, 0])
lam0 = Quaternion([1, 2, 3, 4])
lam0 = Quaternion([-0.235019, -0.144020, 0.502258, 0.819610]) #  GLONASS
lam0 = lam0 * Quaternion([1.0/lam0.getNorm()**0.5, 0, 0, 0])
Nb = 0.35

BASE = 0
    
def lamCircle(phi):
    ans = lam0
    return ans
    
def dlamCircle(phi):
    ans = Quaternion([0, 0, 0, 0])
    return ans
    

def r(phi, ez):
    return 1.0 / (1.0 + ez * cos(phi))

def Nk(phi, k, ez):
    ans = [phi ** k,
           (phi / T) ** k,
           sin((pi * k * phi) / (2 * T)),
           sin((k * phi) / (2 * T)),
           (r(phi, ez) - r(0, ez)) ** k,
           ((r(phi, ez) - r(0, ez))/T) ** k,]
    return ans[BASE]
    
def dNk(phi, k, ez):
    ans = [k * (phi ** (k-1)),
           k * (phi ** (k-1)) / (T ** k),
           ((pi * k) / (2 * T)) * cos((pi * k * phi) / (2 * T)),
           ((k) / (2 * T)) * cos((k * phi) / (2 * T)),
            k * (ez * sin(phi) * (r(phi, ez) ** 2.0)) * ((r(phi, ez) - r(0, ez)) ** (k-1)),
            k * (ez * sin(phi) * (r(phi, ez) ** 2.0)/(T**k)) * ((r(phi, ez) - r(0, ez)) ** (k-1))]
    return ans[BASE]


def lam(a, phi, ez):
    ans = lamCircle(phi)
    for k in range(len(a)):
        # print("ak", a[k])
        # Nkq = Quaternion([Nk(phi, k+1), 0, 0, 0])
        # print("Nkq", Nkq)
        # mul = a[k] * Nkq
        # print("mul", mul)
        ans = ans + a[k] * Quaternion([Nk(phi, k+1, ez), 0, 0, 0])
    return ans

def fvect(phi, ez, s):
    r3 = r(phi, ez) ** 3.0
    omega = Quaternion([0, Nb * r3 * cos(phi), Nb * r3 * sin(phi), 0])
    ans = lamCircle(phi) * omega - Quaternion([2.0, 0, 0, 0]) * dlamCircle(phi)
    
    ans = ans * Quaternion([Nk(phi, s+1, ez), 0, 0, 0])
    return ans
    
def Kvect(phi, ez, s, k):
    # k+1 or k?
    ans = Quaternion([2 * dNk(phi, k+1, ez), 0, 0, 0])
    
    r3 = r(phi, ez) ** 3.0
    omega = Quaternion([0, Nb * r3 * cos(phi), Nb * r3 * sin(phi), 0])
    
    ans = ans - omega * Quaternion([Nk(phi, k+1, ez), 0, 0, 0])
    ans = ans * Quaternion([Nk(phi, s+1, ez), 0, 0, 0])
    return ans
    
# Cauchy
def eqEllipse(lam, phi, ez):
    #print("lam =", lam)
    lam = Quaternion(lam)
    r3 = r(phi, ez) ** 3.0
    omega = Quaternion([0, Nb * r3 * cos(phi), Nb * r3 * sin(phi), 0])
    ans = lam * omega
    ans = ans * Quaternion([0.5, 0, 0, 0])
    return [ans[idx] for idx in range(4)]
    
def genSol(e, M):
    K = [[Quaternion([0, 0, 0, 0]) for _ in range(M)] for _ in range(M)]
    # print("K =", K)
    
    f = [Quaternion([0, 0, 0, 0]) for _ in range(M)]
    # print("f =", f)
    
    for s in range(M):
        phis = (s+1) * T/(M+1)
#         phis = (s+1) * T/M
        # print("phis", phis)
        for k in range(M):
            # pointwise collocation
            K[s][k] = Kvect(phis, e, s, k)
            # print("K[{0}, {1}] = {2}".format(s, k, K[s][k]))
            
        f[s] = [fvect(phis, e, s)]
        # print("f[", s, "] =", f[s])
    # print("f =", f)
    
    aGauss = solveLinear(K, f)
    # for idx, elem in enumerate(aGauss):
    #     print("aGauss[{0}] = {1}".format(idx, list(map(str, elem))))
    a = deepcopy(aGauss[0])
    # print("a =", list(map(str,a)))
    return a

def genError(finish, de):
    res = []
    out = open("res/BASE{0}_big.dat".format(BASE), "w")
    outL = open("res/BASE_Latex{0}_big.dat".format(BASE), "w")
    mx = None
    print("T =", T)
    out.write("T = {0}".format(T))
    # Mmin, Mmax = (2, 9)
    # for m in range (Mmin, Mmax):
    #     out.write(" {0}".format(m))
    # out.write("\n")
    # for m in range (Mmin, Mmax):
    #     outL.write(" {0}".format(m))
    # outL.write("\n")

    for e in np.arange(0, finish+de/2, de):
        print(e)
        out.write("\n{0:.2f}".format(e))
        outL.write("\n${0:.2f}$".format(e))
        for M in range(Mmin, Mmax):
            print("M =", M)
            a = genSol(e, M)
            # print(a)
            phi = np.linspace(0, T, 1001)
            sol = odeint(eqEllipse, [lam0[idx] for idx in range(4)], phi, args = (e,), rtol=1e-15)
            err = []
            for cur in zip(phi, sol):
                l = lam(a, cur[0], e)
                # l = l * Quaternion([1/l.getNorm()**0.5, 0, 0, 0])
                diff = Quaternion(cur[1]) - l
                err.append(diff.getNorm() ** .5)
            mx = max(err)
            print(mx)
            # out.write(" {0:.20f}".format(mx))
            out.write(" {0:.1e}".format(mx))
            latex = " & ${0:.1e}".format(mx) + "}$"
            latex = latex.replace("e-0", " \cdot 10^{-")
            outL.write(latex)
        outL.write("\\\\[0.25em]\n\\hline\n%")
        res.append(mx)
        print()
    out.close()
    outL.close()
    print(res)
    
def genGraph(finish, de, m):
    out = open("res/BASE_Graph{0}_M{1}_big.dat".format(BASE, m), "w")
    # for _ in range(5):
    #     out.write(" {0}".format(0))
    for e in np.arange(0, finish+de/2, de):
        print(e)
        out.write("\n{0:.2f}".format(e))
        print("M =", m)
        a = genSol(e, m)
        # print(a)
        phi = np.linspace(0, T, 1001)
        sol = odeint(eqEllipse, [lam0[idx] for idx in range(4)], phi, args = (e,), rtol=1e-15)
        err = []
        mx, norm = Quaternion([-1, -1, -1, -1]), -1
        for cur in zip(phi, sol):
            l = lam(a, cur[0], e)
            # l = l * Quaternion([1/l.getNorm()**0.5, 0, 0, 0])
            diff = Quaternion(cur[1]) - l
            # print(diff)
            ndiff = diff.getNorm() ** .5
            if ndiff > norm:
                mx = diff
                norm = ndiff
        for idx in range(4):
            out.write(" {0}".format(abs(mx[idx])))
        out.write(" {0:.1e}".format(Quaternion(mx).getNorm() ** .5))
        print(Quaternion(mx).getNorm() ** .5)
    print()
    out.close()
    
    
if __name__ == "__main__":
    print("MVN >>>")    
    # print(r(3.14, EZ))
    # 
    # print(lam([Quaternion([1, 0, 0, 0]), Quaternion([0, 1, 1, 1])], pi/2))
    
        
    
    # lam0[0] = 9
    # print(lam0)
    
    #assert(M == 1)
    print("M =", M)
    
    # a = K[0][0].getInv()
    # print("Inv check (1, 0, 0, 0) =", a * K[0][0])
    # 
    # a = f[0][0] * a 
    # print("a =", a)
    
    a = genSol(EZ, M)
    
    # genError(1e-1, 1e-2)
    
    #genGraph(1e-1, 1e-2, M)
    genGraph(0.51, 0.05, M)
    
    print("lam0 = ", lam0, lam0.getNorm())
    print("lam(0) = ", lam(a, 0, EZ))
    print("lam(T) = ", lam(a, T, EZ), lam(a, T, EZ).getNorm())
    
    #Cauchy
    print("Cauchy >>>")
    phi = np.linspace(0, T, 1001)
    sol = odeint(eqEllipse, [lam0[idx] for idx in range(4)], phi, args = (EZ,), rtol=1e-15)
    # print(sol[:5], "\n\n", sol[-5:])
    
    lastQ = Quaternion(sol[-1])
    print(lastQ, lastQ.getNorm())
    
    print(lam(a, T, EZ), lam(a, T, EZ).getNorm())
    
    err = []
    for cur in zip(phi, sol):
        l = lam(a, cur[0], EZ)
        diff = Quaternion(cur[1]) - l
        err.append(diff.getNorm() ** .5)
        
    # print(err[:5])
    # print(err[-5:])
    
    plt.plot(phi, err, label ='err')
    plt.legend(loc='best')
    plt.xlabel('phi')
    plt.grid()
    plt.show()
    
    print("Cauchy <<<")
    
    print("MVN <<<")    


