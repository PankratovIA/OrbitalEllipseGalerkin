import numpy as np
from MathematicModels import Quaternion
from copy import deepcopy

def makeTriangle(a):
    N = len(a)
    assert(all([len(row) == N+1 for row in a]))
    for r in range(N):
        tmp = a[r][r].getInv()
        a[r] = [x * tmp for x in a[r]]
        for j in range(r+1, N):
            tmp = a[j][r]
            for col in range(r, N+1):
                a[j][col] = a[j][col] - (a[r][col] * tmp)

    return a
    
def getSolution(a):
    """
        a is a triangle
    """
    N = len(a)
    assert(all([len(row) == N+1 for row in a]))
    ans = [a[r][-1] for r in range(N)]
    for r in range(N, -1, -1):
        for prev in range(r+1, N):
            ans[r] = ans[r] - a[r][prev] * ans[prev]
    return ans

def printQuatMatr(a):
    for row in a:
        print("  ".join(map(str,row)))

if __name__ == "__main__":
    print("Gauss >>>")
    a = np.array([[1, 2], [3,4]])
    print("a =", a)
    
    b = np.array([[1], [1]])
    print("b =", b)
    
    x = np.dot(np.linalg.inv(a), b)
    print("x =", x)
    
    aQ = [list(Quaternion([0.2*i+0.8*j+1, 0, 0, 0]) for j in range(2)) for i in range(2)]
    aQ[0].append(Quaternion([1, 0, 0, 0]))
    aQ[1].append(Quaternion([-2, 0, 0, 0]))
    
    print("aQ =")
    printQuatMatr(aQ)
    
#     If the list contains objects and you want to copy them as well,
#     use generic copy.deepcopy()
    aQtr = makeTriangle(deepcopy(aQ))
    print("aQtr =")
    printQuatMatr(aQtr)
    
    x = getSolution(aQtr)
    print("x =", list(map(str, x)))

    npQ = np.array([row[:-1] for row in aQ])
    print("npQ = ")
    printQuatMatr(npQ)
    
    f = np.dot(npQ, np.array([[cur] for cur in x]))
    
    print("f = ")
    printQuatMatr(f)
    
    err = list(map(lambda x: x[0] - x[1], zip(f, [deepcopy(aQ[r][-1]) for r in range(len(aQ))])))
    print("err = ")
    printQuatMatr(err)
    print("Gauss <<<") 