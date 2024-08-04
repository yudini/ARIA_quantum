import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import ResourceCounter, ClassicalSimulator
from projectq.meta import Compute, Uncompute, Dagger
def CNOT8bit(a, an, b, bn):
    for i in range(8):
        CNOT | (a[127 - an * 8 - i], b[127 - bn * 8 - i])

def DL(eng,a):
    y = eng.allocate_qureg(128)

    CNOT8bit(a, 3, y, 0)
    CNOT8bit(a, 2, y, 1)
    CNOT8bit(a, 1, y, 2)
    CNOT8bit(a, 0, y, 3)
    CNOT8bit(a, 5, y, 4)
    CNOT8bit(a, 4, y, 5)
    CNOT8bit(a, 7, y, 6)
    CNOT8bit(a, 6, y, 7)

    CNOT8bit(a, 15, y, 8)
    CNOT8bit(a, 14, y, 9)
    CNOT8bit(a, 13, y, 10)
    CNOT8bit(a, 12, y, 11)
    CNOT8bit(a, 11, y, 12)
    CNOT8bit(a, 10, y, 13)
    CNOT8bit(a, 9, y, 14)
    CNOT8bit(a, 8, y, 15)

#######################################

    CNOT8bit(a, 4, y, 0)
    CNOT8bit(a, 5, y, 1)
    CNOT8bit(a, 6, y, 2)
    CNOT8bit(a, 7, y, 3)
    CNOT8bit(a, 0, y, 4)
    CNOT8bit(a, 1, y, 5)
    CNOT8bit(a, 2, y, 6)
    CNOT8bit(a, 3, y, 7)

    CNOT8bit(a, 13, y, 8)
    CNOT8bit(a, 12, y, 9)
    CNOT8bit(a, 15, y, 10)
    CNOT8bit(a, 14, y, 11)
    CNOT8bit(a, 9, y, 12)
    CNOT8bit(a, 8, y, 13)
    CNOT8bit(a, 11, y, 14)
    CNOT8bit(a, 10, y, 15)

    #######################################

    CNOT8bit(a, 6, y, 0)
    CNOT8bit(a, 7, y, 1)
    CNOT8bit(a, 4, y, 2)
    CNOT8bit(a, 5, y, 3)
    CNOT8bit(a, 2, y, 4)
    CNOT8bit(a, 3, y, 5)
    CNOT8bit(a, 0, y, 6)
    CNOT8bit(a, 1, y, 7)

    CNOT8bit(a, 10, y, 8)
    CNOT8bit(a, 11, y, 9)
    CNOT8bit(a, 8, y, 10)
    CNOT8bit(a, 9, y, 11)
    CNOT8bit(a, 12, y, 12)
    CNOT8bit(a, 13, y, 13)
    CNOT8bit(a, 14, y, 14)
    CNOT8bit(a, 15, y, 15)

    #######################################

    CNOT8bit(a, 14, y, 0)
    CNOT8bit(a, 15, y, 1)
    CNOT8bit(a, 12, y, 2)
    CNOT8bit(a, 13, y, 3)
    CNOT8bit(a, 11, y, 4)
    CNOT8bit(a, 10, y, 5)
    CNOT8bit(a, 9, y, 6)
    CNOT8bit(a, 8, y, 7)

    CNOT8bit(a, 7, y, 8)
    CNOT8bit(a, 6, y, 9)
    CNOT8bit(a, 5, y, 10)
    CNOT8bit(a, 4, y, 11)
    CNOT8bit(a, 2, y, 12)
    CNOT8bit(a, 3, y, 13)
    CNOT8bit(a, 0, y, 14)
    CNOT8bit(a, 1, y, 15)

    #######################################

    CNOT8bit(a, 13, y, 0)
    CNOT8bit(a, 12, y, 1)
    CNOT8bit(a, 15, y, 2)
    CNOT8bit(a, 14, y, 3)
    CNOT8bit(a, 8, y, 4)
    CNOT8bit(a, 9, y, 5)
    CNOT8bit(a, 10, y, 6)
    CNOT8bit(a, 11, y, 7)

    CNOT8bit(a, 4, y, 8)
    CNOT8bit(a, 5, y, 9)
    CNOT8bit(a, 6, y, 10)
    CNOT8bit(a, 7, y, 11)
    CNOT8bit(a, 1, y, 12)
    CNOT8bit(a, 0, y, 13)
    CNOT8bit(a, 3, y, 14)
    CNOT8bit(a, 2, y, 15)

    #######################################

    CNOT8bit(a, 9, y, 0)
    CNOT8bit(a, 8, y, 1)
    CNOT8bit(a, 11, y, 2)
    CNOT8bit(a, 10, y, 3)
    CNOT8bit(a, 15, y, 4)
    CNOT8bit(a, 14, y, 5)
    CNOT8bit(a, 13, y, 6)
    CNOT8bit(a, 12, y, 7)

    CNOT8bit(a, 1, y, 8)
    CNOT8bit(a, 0, y, 9)
    CNOT8bit(a, 3, y, 10)
    CNOT8bit(a, 2, y, 11)
    CNOT8bit(a, 7, y, 12)
    CNOT8bit(a, 6, y, 13)
    CNOT8bit(a, 5, y, 14)
    CNOT8bit(a, 4, y, 15)

    #######################################

    CNOT8bit(a, 8, y, 0)
    CNOT8bit(a, 9, y, 1)
    CNOT8bit(a, 10, y, 2)
    CNOT8bit(a, 11, y, 3)
    CNOT8bit(a, 14, y, 4)
    CNOT8bit(a, 15, y, 5)
    CNOT8bit(a, 12, y, 6)
    CNOT8bit(a, 13, y, 7)

    CNOT8bit(a, 0, y, 8)
    CNOT8bit(a, 1, y, 9)
    CNOT8bit(a, 2, y, 10)
    CNOT8bit(a, 3, y, 11)
    CNOT8bit(a, 6, y, 12)
    CNOT8bit(a, 7, y, 13)
    CNOT8bit(a, 4, y, 14)
    CNOT8bit(a, 5, y, 15)


    return y

def Diffusion_old(eng,a):
    L = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,
         1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
         1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
         0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0,
         1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0,
         1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0,
         0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1]

    U = [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0,
         0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1,
         0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1,
         0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1,
         0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1,
         0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

    n = 16
    for i in range(n - 1):
        for j in range(n - 1 - i):
            if (U[(i * n) + 1 + i + j] == 1):
                CNOT8bit(a, 1 + i + j, a, i)

    for i in range(n - 1):
        for j in range(n - 1, 0 + i, -1):
            if (L[n * (n - 1 - i) + (n - 1 - j)] == 1):
                CNOT8bit(a, n - 1 - j, a, n - 1 - i)

    out = []
    for i in range(4):
        for j in range(8):
            out.append(a[8 * i + j])
    for j in range(8):
        out.append(a[8 * 6 + j])
    for j in range(8):
        out.append(a[8 * 7 + j])
    for j in range(8):
        out.append(a[8 * 4 + j])
    for j in range(8):
        out.append(a[8 * 5 + j])
    for j in range(8):
        out.append(a[8 * 9 + j])
    for j in range(8):
        out.append(a[8 * 8 + j])
    for j in range(8):
        out.append(a[8 * 11 + j])
    for j in range(8):
        out.append(a[8 * 10 + j])
    for i in range(15, 11, -1):
        for j in range(8):
            out.append(a[8 * i + j])
    return out



def Diffusion(eng, a):
    CNOT8bit(a, 3, a, 9)
    CNOT8bit(a, 13, a, 3)
    CNOT8bit(a, 6, a, 8)
    CNOT8bit(a, 4, a, 14)
    CNOT8bit(a, 12, a, 6)
    CNOT8bit(a, 9, a, 13)
    CNOT8bit(a, 1, a, 11)
    CNOT8bit(a, 9, a, 1)
    CNOT8bit(a, 8, a, 12)
    CNOT8bit(a, 6, a, 4)
    CNOT8bit(a, 7, a, 13)
    CNOT8bit(a, 13, a, 6)
    CNOT8bit(a, 14, a, 7)
    CNOT8bit(a, 8, a, 7)
    CNOT8bit(a, 10, a, 14)
    CNOT8bit(a, 4, a, 10)
    CNOT8bit(a, 2, a, 12)
    CNOT8bit(a, 1, a, 4)
    CNOT8bit(a, 15, a, 1)
    CNOT8bit(a, 5, a, 15)
    CNOT8bit(a, 10, a, 5)
    CNOT8bit(a, 6, a, 10)
    CNOT8bit(a, 3, a, 2)
    CNOT8bit(a, 0, a, 14)
    CNOT8bit(a, 2, a, 8)
    CNOT8bit(a, 11, a, 2)
    CNOT8bit(a, 12, a, 9)
    CNOT8bit(a, 2, a, 12)
    CNOT8bit(a, 15, a, 8)
    CNOT8bit(a, 1, a, 0)
    CNOT8bit(a, 7, a, 3)
    CNOT8bit(a, 9, a, 1)
    CNOT8bit(a, 7, a, 13)
    CNOT8bit(a, 9, a, 7)
    CNOT8bit(a, 5, a, 11)
    CNOT8bit(a, 11, a, 15)
    CNOT8bit(a, 14, a, 3)
    CNOT8bit(a, 14, a, 11)
    CNOT8bit(a, 1, a, 5)
    CNOT8bit(a, 11, a, 4)
    CNOT8bit(a, 11, a, 1)
    CNOT8bit(a, 3, a, 9)
    CNOT8bit(a, 0, a, 14)
    CNOT8bit(a, 10, a, 0)
    CNOT8bit(a, 4, a, 10)
    CNOT8bit(a, 6, a, 2)
    CNOT8bit(a, 8, a, 6)

    out = []
    for j in range(8):
        out.append(a[8 * (15 - 5) + j])

    for j in range(8):
        out.append(a[8 * (15 - 4) + j])

    for j in range(8):
        out.append(a[8 * (15 - 3) + j])

    for j in range(8):
        out.append(a[8 * (15 - 2) + j])

    for j in range(8):
        out.append(a[8 * (15 - 7) + j])

    for j in range(8):
        out.append(a[8 * (15 - 8) + j])

    for j in range(8):
        out.append(a[8 * (15 - 11) + j])

    for j in range(8):
        out.append(a[8 * (15 - 0) + j])

    for j in range(8):
        out.append(a[8 * (15 - 12) + j])

    for j in range(8):
        out.append(a[8 * (15 - 9) + j])

    for j in range(8):
        out.append(a[8 * (15 - 14) + j])

    for j in range(8):
        out.append(a[8 * (15 - 1) + j])

    for j in range(8):
        out.append(a[8 * (15 - 10) + j])

    for j in range(8):
        out.append(a[8 * (15 - 15) + j])

    for j in range(8):
        out.append(a[8 * (15 - 6) + j])

    for j in range(8):
        out.append(a[8 * (15 - 13) + j])
    return out

def print_state(eng, b, n):  # n = /4
    All(Measure) | b
    print('0x', end='')
    print_hex(eng, b, n)
    print('\n')


def print_hex(eng, qubits, n):
    for i in reversed(range(n)):
        temp = 0
        temp = temp + int(qubits[4 * i + 3]) * 8
        temp = temp + int(qubits[4 * i + 2]) * 4
        temp = temp + int(qubits[4 * i + 1]) * 2
        temp = temp + int(qubits[4 * i])

        temp = hex(temp)
        y = temp.replace("0x", "")
        print(y, end='')

def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]


def ARIA128(eng):
    r = 12  # round number
    bit = 8
    n = 128

    PT = eng.allocate_qureg(n)

    if (resource_check != 1):
        Round_constant_XOR(eng, PT, P, 128)

    #PT = Diffusion(eng, PT)
    #PT = DL(eng, PT)
    PT = Diffusion_old(eng,PT)

    if(resource_check!=1):
        print_state(eng,PT,32)


MK = 0x00112233445566778899aabbccddeeff
# P = 0x11111111aaaaaaaa11111111bbbbbbbb
# P = 0x11111111aaaaaaaa11111111abcdefff
P = 0x517cc1b727220a94fe13abe8fa9a6ee0
NCT = 1
resource_check = 0
AND_check = 0

# 결과값 확인
classic = ClassicalSimulator()
eng = MainEngine(classic)
print("#### Ciphertext ####")
ARIA128(eng)
eng.flush()
print('\n')

# 분해 전
NCT = 1
resource_check = 1
Resource = ResourceCounter()
eng = MainEngine(Resource)
print("#### Quantum Resources costs (NCT Level) ####")
ARIA128(eng)
print(Resource)
eng.flush()
print('\n')

# 분해 후
NCT = 0
Resource = ResourceCounter()
eng = MainEngine(Resource)
print("#### Decomposed Quantum Resources costs ####")
ARIA128(eng)
print(Resource)
print('\n')