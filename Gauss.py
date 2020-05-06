import numpy as np
from MathematicModels import Quaternion

def makeTriangle(a):
    N = len(a)
    assert(all([len(row) == N+1 for row in a]))
    for r in range(N-1):
        tmp = a[r][r].getInv()
        #print("a[{0}] =".format(row))
        #printQuatMatr(a[row])
        a[r] = [x * tmp for x in a[r]]
        #print("a[{0}] =".format(row))
        #printQuatMatr(a[row])
        for j in range(r+1, N):
            tmp = a[j][r]
#             print("tmp =", str(tmp))
            for col in range(r, N+1):
#                 print("a[j][col] a[r][col] tmp", str(a[j][col]), str(a[r][col]), str(tmp))
#                 print("a[r][col] * tmp", str(a[r][col] * tmp))
                z = a[j][col] - (a[r][col] * tmp)
#                 print("z", str(z))
#                 print(" a[j][col] >>>", str(a[j][col]))
                a[j][col] = z#Quaternion([z[0], z[1], z[2], z[3]])
#                 print(" a[j][col] <<<", str(a[j][col]))
                
#     for i in range(N):
#         for j in range(N+1):
#             print("i = {0}, j = {1} {2}".format(i, j, str(a[i][j])))
#     printQuatMatr(a)
                
    return list(a)

def printQuatMatr(a):
    for row in a:
        #print("row =", row)
        print("  ".join(map(str,row)))

if __name__ == "__main__":
    print("Gauss >>>")
    a = np.array([[1, 2], [3,4]])
    print("a =", a)
    
    b = np.array([[1], [1]])
    print("b =", b)
    
    x = np.dot(np.linalg.inv(a), b)
    print("x =", x)
    
    aQ = [list(Quaternion([i+j+1, 0, 0, 0]) for j in range(2)) for i in range(2)]
    aQ[0].append(Quaternion([1, 0, 0, 0]))
    aQ[1].append(Quaternion([1, 0, 0, 0]))
#     print("aQ =", aQ)
#     for i in range(2):
#         print("---")
#         for j in range(3):
#             print("i = {0}, j = {1}, {2}".format(i, j, aQ[i][j]))
    
    print("aQ =")
    printQuatMatr(aQ)
    
    aQtr = makeTriangle(aQ)
    print("aQtr =")
    printQuatMatr(aQtr)
    

    print("Gauss <<<") 