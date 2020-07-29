from math import sin, cos, pi
import numpy as np
from MathematicModels import Quaternion
from Gauss import solveLinear
from copy import deepcopy
from scipy.integrate import quad, odeint, ode
import matplotlib.pyplot as plt

EZ = 0.01
T = pi/2
M = 5
# lam0 = Quaternion([1, 0, 0, 0])
# lam0 = Quaternion([0.8, 0, 0.6, 0])
lam0 = Quaternion([1, 2, 3, 4])
lam0 = Quaternion([-0.235019, -0.144020, 0.502258, 0.819610]) #  GLONASS
lam0 = lam0 * Quaternion([1.0/lam0.getNorm()**0.5, 0, 0, 0])
Nb = 0.35

BASE = 2
    
def lamCircle(phi):
    omega = Quaternion([0, Nb, 0, 1.0])
    om = (1 + Nb * Nb) ** 0.5
    
    omega = omega * Quaternion([sin(0.5 * om * phi)/om, 0, 0, 0])
    ans = Quaternion([cos(0.5 * om * phi), 0, 0, 0])
    ans = ans + omega
    ans = lam0 * ans
    # ans = lam0
    return ans
    
def dlamCircle(phi):
    omega = Quaternion([0, Nb, 0, 1.0])
    om = (1 + Nb * Nb) ** 0.5
    omega = omega * Quaternion([0.5*cos(0.5 * om * phi), 0, 0, 0])
    ans = Quaternion([-0.5 * om * sin(0.5 * om * phi), 0, 0, 0])
    ans = ans + omega
    ans = lam0 * ans
    # ans = Quaternion([0, 0, 0, 0])
    return ans
    

def r(phi, ez):
    return 1.0 / (1.0 + ez * cos(phi))

def Nk(phi, k):
    ans = [phi ** k,
           (phi / T) ** k,
           sin((pi * k * phi) / (2 * T)),
           sin((k * phi) / (2 * T))]#, (r(phi) - r(0)) ** k]
    return ans[BASE]
    
def dNk(phi, k):
    ans = [k * (phi ** (k-1)),
           k * (phi ** (k-1)) / (T ** k),
           ((pi * k) / (2 * T)) * cos((pi * k * phi) / (2 * T)),
           ((k) / (2 * T)) * cos((k * phi) / (2 * T))]
           # k * (ez * sin(phi) * (r(phi) ** 2.0)) * ((r(phi) - r(0)) ** (k-1))
    return ans[BASE]


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

def fvect(phi, ez, s):
    r3 = r(phi, ez) ** 3.0
    omega = Quaternion([0, Nb * r3, 0, 1.0])
    ans = lamCircle(phi) * omega - Quaternion([2.0, 0, 0, 0]) * dlamCircle(phi)
    
    ans = ans * Quaternion([Nk(phi, s+1), 0, 0, 0])
    # short formula from paper
    short = lamCircle(phi) * Quaternion([0, Nb*(r3 - 1.0), 0, 0])
    short = short * Quaternion([Nk(phi, s+1), 0, 0, 0])
    # print("ans -- short", ans[idx], short[idx], abs(ans[idx] - short[idx]))
    return short
    
def Kvect(phi, ez, s, k):
    # k+1 or k?
    ans = Quaternion([2 * dNk(phi, k+1), 0, 0, 0])
    
    r3 = r(phi, ez) ** 3.0
    omega = Quaternion([0, Nb * r3, 0, 1.0])
    
    ans = ans - omega * Quaternion([Nk(phi, k+1), 0, 0, 0])
    ans = ans * Quaternion([Nk(phi, s+1), 0, 0, 0])
    return ans
    
# Cauchy
def eqEllipse(lam, phi, ez):
    #print("lam =", lam)
    lam = Quaternion(lam)
    r3 = r(phi, ez) ** 3.0
    omega = Quaternion([0, Nb * r3, 0, 1.0])
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
    out = open("res/BASE{0}.dat".format(BASE, T), "w")
    outL = open("res/BASE_Latex{0}.dat".format(BASE, T), "w")
    mx = None
    print("T =", T)
    out.write("T = {0}".format(T))
    Mmin, mMax = (8, 13)
    for m in range (Mmin, mMax):
        out.write(" {0}".format(m))
    out.write("\n")
    for m in range (Mmin, mMax):
        outL.write(" {0}".format(m))
    outL.write("\n")

    for e in np.arange(de, finish+de/2, de):
        print(e)
        out.write("\n{0:.2f}".format(e))
        outL.write("\n${0:.2f}$".format(e))
        for M in range(8, 13):
            print("M =", M)
            a = genSol(e, M)
            # print(a)
            phi = np.linspace(0, T, 1001)
            sol = odeint(eqEllipse, [lam0[idx] for idx in range(4)], phi, args = (e,), rtol=1e-15)
            err = []
            for cur in zip(phi, sol):
                l = lam(a, cur[0])
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
    
    genError(1e-1, 1e-2)
    
    print("lam0 = ", lam0, lam0.getNorm())
    print("lam(0) = ", lam(a, 0))
    print("lam(T) = ", lam(a, T), lam(a, T).getNorm())
    
    #Cauchy
    print("Cauchy >>>")
    phi = np.linspace(0, T, 1001)
    sol = odeint(eqEllipse, [lam0[idx] for idx in range(4)], phi, args = (EZ,), rtol=1e-15)
    # print(sol[:5], "\n\n", sol[-5:])
    
    lastQ = Quaternion(sol[-1])
    print(lastQ, lastQ.getNorm())
    
    print(lam(a, T), lam(a, T).getNorm())
    
    err = []
    for cur in zip(phi, sol):
        l = lam(a, cur[0])
        diff = Quaternion(cur[1]) - l
        err.append(diff.getNorm() ** .5)
        
    # print(err[:5])
    # print(err[-5:])
    
    # plt.plot(phi, err, label ='err')
    # plt.legend(loc='best')
    # plt.xlabel('phi')
    # plt.grid()
    # plt.show()
    
    print("Cauchy <<<")
    
    print("MVN <<<")    


