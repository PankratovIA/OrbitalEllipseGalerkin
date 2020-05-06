import numpy as np
from MathematicModels import Quaternion

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
    
    aQ = [list(Quaternion([i+j+2, 0, 0, 0]) for j in range(2)) for i in range(2)]
    aQ[0].append(Quaternion([1, 0, 0, 0]))
    aQ[1].append(Quaternion([2, 0, 0, 0]))
    
    print("aQ =")
    printQuatMatr(aQ)
    
    aQtr = makeTriangle(aQ)
    print("aQtr =")
    printQuatMatr(aQtr)
    
    x = getSolution(aQtr)
    print("x =", list(map(str, x)))
#     printQuatMatr(x)
    

    print("Gauss <<<") 