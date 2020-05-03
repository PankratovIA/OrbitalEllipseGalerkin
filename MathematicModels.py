#TuraevA - 411
import math as mt
import traceback


class Quaternion:

    def __init__(self, variable_list):
        checkList = False
        if isinstance(variable_list, type(list())):
            for i in variable_list:
                if not (isinstance(i, int) or isinstance(i, float)):
                    checkList = True
                    break
        if len(variable_list) != 4 or checkList:
            raise SystemExit("Wrong Quaternion", variable_list)
        self.var = variable_list

    def getConjugate(self):
        return Quaternion([self.var[0], -self.var[1], -self.var[2], -self.var[3]])
        
    def getInv(self):
        conj = self.getConjugate()
        # norm = self.getNorm()
        return conj * Quaternion([self.getNorm(), 0, 0, 0])

    def __add__(self, other):
        # print("other", other)
        return Quaternion([self.var[i] + other.var[i] for i in range(4)])

    def __iadd__(self, other):
        for i in range(4):
            self.var[i] += other.var[i]

    def __sub__(self, other):
        return Quaternion([self.var[i] - other.var[i] for i in range(4)])

    def __isub__(self, other):
        for i in range(4):
            self.var[i] -= other.var[i]

    def __mul__(self, other):
        scalar_value = [self.var[0] * other.var[0], 0, 0, 0]
        for i in range(1, 4):
            scalar_value[0] -= self.var[i] * other.var[i]
        mul1 = [self.var[0] * other.var[i] for i in range(4)]
        mul2 = [other.var[0] * self.var[i] for i in range(4)]
        mul1[0] = mul2[0] = 0
        vecMul = [0, self.var[2] * other.var[3] - self.var[3] * other.var[2],
                  -(self.var[1] * other.var[3] - self.var[3] * other.var[1]),
                  self.var[1] * other.var[2] - self.var[2] * other.var[1]]
        return Quaternion([scalar_value[i] + mul1[i] + mul2[i] + vecMul[i] for i in range(4)])

    def __truediv__(self, other):
        try:
            result = self * Quaternion(other.getConjugate())
            norm = other.getNorm()
            return Quaternion([result[i] / norm for i in range(4)])
        except ZeroDivisionError as exc:
            raise SystemExit(traceback.print_tb(exc.__traceback__))

    def getTrigonometrical(self):
        temp1 = 0
        temp2 = 0
        for i in range(4):
            temp1 += self.var[i] ** 2
            if i > 0:
                temp2 += self.var[i] ** 2
        fullRoot = mt.sqrt(temp1)
        root = mt.sqrt(temp2)
        # angleCos = mt.acos(self.var[0] / fullRoot)
        # angleSin = mt.asin(root / fullRoot)
        angleTang = mt.atan(root / self.var[0])
        ksi = [self.var[i] / root for i in range(1, 4)]
        return '{0:.5}(cos({1:.5}) + ({2[0]:.5} i1 + {2[1]:.5} i2 + {2[2]:.5} i3)sin({3:.5}))'. \
            format(fullRoot, mt.degrees(angleTang), ksi, mt.degrees(angleTang))

    def __getitem__(self, idx):
        assert(0<=idx<4)
        return self.var[idx]
        
    def __setitem__(self, idx, elem):
        assert(0<=idx<4)
        self.var[idx] = elem

    def getNorm(self):
        return (self * self.getConjugate())[0]

    def getVar(self):
        return self.var

    def __str__(self):
        return '{0[0]} + {0[1]}*i1 + {0[2]}*i2 + {0[3]}*i3'.format(self.var)
    
        
    

def turn_coordinate(angle: int, r=(0, 0, 0), e=(0, 0, 0)):
    r = Quaternion([0] + r)
    angle = mt.radians(angle)
    lambd = Quaternion([mt.cos(angle / 2), [var * mt.sin(angle / 2) for var in e]])
    new_r = [0, 0, 0]
    quatern = Quaternion(lambd.getConjugate()) * Quaternion(r * lambd)
    for i in range(1, 4):
        new_r[i - 1] = round(quatern[i], 10)
    return new_r


def turn_vector(angle: int, r=(0, 0, 0), e=(0, 0, 0)):
    r = Quaternion([0] + r)
    angle = mt.radians(angle)
    lambd = Quaternion([mt.cos(angle / 2), [var * mt.sin(angle / 2) for var in e]])
    new_r = [0, 0, 0]
    quatern = Quaternion(lambd * r) * Quaternion(lambd.getConjugate())
    for i in range(1, 4):
        new_r[i - 1] = round(quatern[i], 7)
    return new_r


def main():
    a = Quaternion([-2, 0, 3, 1])
    b = Quaternion([0, 2, 1, -3])
    # a = Quaternion([0.5, 0, 0, mt.sqrt(3) / 2])
    # b = Quaternion([mt.sqrt(3) / 2, 0, -0.5, 0])

    print("a =", str(a))
    print("b =", str(b))
    print("difference: a - b =", a - b)
    print("addition: a + b =", a + b)
    print("multiplication: a * b =", a * b)
    # print("check multiplication with norm: a * b =", Quaternion(a * b).getNorm())
    # print("division: a / b =", a / b)
    # print("check division with norm: a / b =", Quaternion(a / b).getNorm())
    # print("check division: a / b * b =", Quaternion(a / b) * b)
    print("trigonometrical form a: ", a.getTrigonometrical())
    print("a conjure = ", a.getConjugate())
    print("b conjure = ", b.getConjugate())
    print("norm a =", a.getNorm())
    print("norm b =", b.getNorm())
    var_turn = turn_coordinate(90, [0, 0, 1], [0, 1, 0])
    print('turn coordinate system 1 = ', var_turn)
    print('turn coordinate system 2 = ', turn_coordinate(90, var_turn, [0, 0, 1]))
    var_turn = turn_vector(90, [0, 0, 1], [0, 1, 0])
    print('turn vector 1 = ', var_turn)
    print('turn vector 2 = ', turn_vector(90, var_turn, [0, 0, 1]))


if __name__ == '__main__':
    main()
