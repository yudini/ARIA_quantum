import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import ResourceCounter, ClassicalSimulator
from projectq.meta import Compute, Uncompute, Dagger


def main(eng):
    a= eng.allocate_qureg(8)
    ancilla = eng.allocate_qureg(38)
    #s = eng.allocate_qureg(8)
    vect_a = 0x63
    vect_b = 0xe2
    P = 0xa9

    flag =0
    round =1

    if (resource_check != 1):
        Round_constant_XOR(eng, a, P, 8)

    #a= Sbox1(eng,a,vect_a,8,ancilla)
    #a= Sbox1Inv(eng,a,vect_a,8,ancilla)
    # a= Sbox2(eng,a,vect_b,8,ancilla)
    a= Sbox2Inv(eng,a,vect_b,8,ancilla)


    if(resource_check !=1):
        print_state(eng,a,2)

def Sbox1(eng, a, vect_a, n, ancilla):
    a1 = eng.allocate_qureg(n)

    res = getInverse(eng, a, a1, ancilla, n)
    res = MatrixProductS1(eng, res, vect_a, n)

    return res

def Sbox1Inv(eng, a, vect_a, n, ancilla):
    a1 = eng.allocate_qureg(n)

    a = MatrixProductS1Inv(eng, a, vect_a, n)
    res = getInverse(eng, a, a1, ancilla, n)

    return res

def Sbox2(eng, a, vect_b, n, ancilla):
    a1 = eng.allocate_qureg(n)

    res = getInverse(eng, a, a1, ancilla, n)
    res = MatrixProductS2(eng, res, vect_b, n)

    return res


def Sbox2Inv(eng, a, vect_b, n, ancilla):
    a1 = eng.allocate_qureg(n)

    a = MatrixProductS2Inv(eng, a, vect_b, n)
    res = getInverse(eng, a, a1, ancilla, n)

    return res

# yang
def getInverse(eng, a, a1, ancilla, n):  # n=8
    count = 0

    copy(eng, a, a1, n)
    a1 = Squaring(eng, a1, n)  # a1 = a^2

    a2 = []
    a2, count, ancilla = recursive_karatsuba(eng, a, a1, n, count, ancilla)
    a2 = Reduction(eng, a2)

    a1 = Squaring(eng, a1, n)
    a1 = Squaring(eng, a1, n)
    a1 = Squaring(eng, a1, n)
    a1 = Squaring(eng, a1, n)
    a1 = Squaring(eng, a1, n)  # ^64

    count = 0
    a3 = []
    a3, count, ancilla = recursive_karatsuba(eng, a2, a1, n, count, ancilla)
    a3 = Reduction(eng, a3)

    a2 = Squaring(eng, a2, n)
    a2 = Squaring(eng, a2, n)

    count = 0
    a4 = []
    a4, count, ancilla = recursive_karatsuba(eng, a2, a3, n, count, ancilla)
    a4 = Reduction(eng, a4)

    a2 = Squaring(eng, a2, n)
    a2 = Squaring(eng, a2, n)

    count = 0
    a5 = []
    a5, count, ancilla = recursive_karatsuba(eng, a2, a4, n, count, ancilla)

    a5 = Reduction(eng, a5)
    a5 = Squaring(eng, a5, n)

    return a5

def MatrixProductS1(eng, inv, b, n):
    out = eng.allocate_qureg(8)

    CNOT | (inv[0], out[0])
    CNOT | (inv[1], out[1])
    CNOT | (inv[2], out[2])
    CNOT | (inv[3], out[3])
    CNOT | (inv[4], out[4])
    CNOT | (inv[5], out[5])
    CNOT | (inv[6], out[6])
    CNOT | (inv[7], out[7])

    CNOT | (inv[7], out[0])
    CNOT | (inv[0], out[1])
    CNOT | (inv[1], out[2])
    CNOT | (inv[2], out[3])
    CNOT | (inv[3], out[4])
    CNOT | (inv[4], out[5])
    CNOT | (inv[5], out[6])
    CNOT | (inv[6], out[7])

    CNOT | (inv[6], out[0])
    CNOT | (inv[7], out[1])
    CNOT | (inv[0], out[2])
    CNOT | (inv[1], out[3])
    CNOT | (inv[2], out[4])
    CNOT | (inv[3], out[5])
    CNOT | (inv[4], out[6])
    CNOT | (inv[5], out[7])

    CNOT | (inv[5], out[0])
    CNOT | (inv[6], out[1])
    CNOT | (inv[7], out[2])
    CNOT | (inv[0], out[3])
    CNOT | (inv[1], out[4])
    CNOT | (inv[2], out[5])
    CNOT | (inv[3], out[6])
    CNOT | (inv[4], out[7])


    CNOT | (inv[4], out[0])
    CNOT | (inv[5], out[1])
    CNOT | (inv[6], out[2])
    CNOT | (inv[7], out[3])
    CNOT | (inv[0], out[4])
    CNOT | (inv[1], out[5])
    CNOT | (inv[2], out[6])
    CNOT | (inv[3], out[7])

    Round_constant_XOR(eng, out, b, n) # out = out ^ b

    return out

def MatrixProductS1Inv(eng, a, b, n):
    Round_constant_XOR(eng, a, b, n) # a = a ^ b

    out = eng.allocate_qureg(8)

    CNOT | (a[2], out[0])
    CNOT | (a[3], out[1])
    CNOT | (a[4], out[2])
    CNOT | (a[5], out[3])
    CNOT | (a[6], out[4])
    CNOT | (a[7], out[5])
    CNOT | (a[0], out[6])
    CNOT | (a[1], out[7])

    CNOT | (a[5], out[0])
    CNOT | (a[6], out[1])
    CNOT | (a[7], out[2])
    CNOT | (a[0], out[3])
    CNOT | (a[1], out[4])
    CNOT | (a[2], out[5])
    CNOT | (a[3], out[6])
    CNOT | (a[4], out[7])


    CNOT | (a[7], out[0])
    CNOT | (a[0], out[1])
    CNOT | (a[1], out[2])
    CNOT | (a[2], out[3])
    CNOT | (a[3], out[4])
    CNOT | (a[4], out[5])
    CNOT | (a[5], out[6])
    CNOT | (a[6], out[7])


    return out

def MatrixProductS2(eng, inv, b, n):
    out = eng.allocate_qureg(8)

    CNOT | (inv[0], out[7])
    CNOT | (inv[1], out[2])
    CNOT | (inv[2], out[1])
    CNOT | (inv[4], out[5])
    CNOT | (inv[6], out[3])
    CNOT | (inv[7], out[0])

    #######################################

    CNOT | (inv[1], out[7])
    CNOT | (inv[2], out[2])
    CNOT | (inv[3], out[1])
    CNOT | (inv[5], out[5])
    CNOT | (inv[6], out[6])
    CNOT | (inv[7], out[3])

    #######################################

    CNOT | (inv[0], out[2])
    CNOT | (inv[1], out[3])
    CNOT | (inv[2], out[7])
    CNOT | (inv[4], out[1])
    CNOT | (inv[5], out[0])
    CNOT | (inv[6], out[5])
    CNOT | (inv[7], out[4])

    #######################################

    CNOT | (inv[1], out[5])
    CNOT | (inv[3], out[7])
    CNOT | (inv[4], out[2])
    CNOT | (inv[5], out[1])
    CNOT | (inv[6], out[0])
    CNOT | (inv[7], out[6])

    #######################################

    CNOT | (inv[0], out[5])
    CNOT | (inv[1], out[6])
    CNOT | (inv[3], out[0])
    CNOT | (inv[5], out[7])
    CNOT | (inv[6], out[1])
    CNOT | (inv[7], out[2])

    #######################################

    CNOT | (inv[0], out[3])
    CNOT | (inv[1], out[0])
    CNOT | (inv[2], out[6])
    CNOT | (inv[5], out[2])
    CNOT | (inv[6], out[7])
    CNOT | (inv[7], out[1])

    #######################################

    CNOT | (inv[1], out[4])
    CNOT | (inv[6], out[4])


    Round_constant_XOR(eng, out, b, n)  # out = out ^ b

    return out

def MatrixProductS2Inv(eng, a, b, n):
    Round_constant_XOR(eng, a, b, n)  # a = a ^ b

    out = eng.allocate_qureg(8)


    CNOT | (a[0], out[6])
    CNOT | (a[1], out[4])
    CNOT | (a[2], out[5])
    CNOT | (a[3], out[0])

    #######################################

    CNOT | (a[0], out[4])
    CNOT | (a[1], out[3])
    CNOT | (a[2], out[6])
    CNOT | (a[4], out[0])
    CNOT | (a[6], out[2])

    #######################################

    CNOT | (a[0], out[3])
    CNOT | (a[2], out[4])
    CNOT | (a[3], out[6])
    CNOT | (a[4], out[5])
    CNOT | (a[6], out[1])
    CNOT | (a[7], out[7])

    #######################################

    CNOT | (a[1], out[5])
    CNOT | (a[4], out[6])
    CNOT | (a[5], out[1])
    CNOT | (a[6], out[7])
    CNOT | (a[7], out[3])

    #######################################

    CNOT | (a[0], out[7])
    CNOT | (a[2], out[1])
    CNOT | (a[4], out[4])
    CNOT | (a[5], out[6])
    CNOT | (a[6], out[3])
    CNOT | (a[7], out[5])

    #######################################

    CNOT | (a[2], out[3])
    CNOT | (a[3], out[7])
    CNOT | (a[4], out[2])
    CNOT | (a[5], out[4])
    CNOT | (a[6], out[5])
    CNOT | (a[7], out[6])

    return out

def recursive_karatsuba(eng, a, b, n, count, ancilla): #n=4

    if(n==1):
        c = eng.allocate_qubit()
        Toffoli_gate(eng, a, b, c)

        return c, count, ancilla

    c_len = 3**math.log(n, 2) #9 #3
    r_low = n//2    #2 #1

    if(n%2!=0):
        r_low = r_low +1 # n=3 -> 2, n=4 -> 2

    r_a = []
    r_b = []

    # Provide rooms and prepare operands
    r_a = ancilla[count:count + r_low]
    #print(count, count + r_low)

    #r_a = room(eng, r_low) #2qubits for r

    count = count + r_low

    r_b = ancilla[count:count + r_low]
    #print(count, count + r_low)

    #r_b = room(eng, r_low) #2qubits for r

    count = count + r_low

    with Compute(eng):
        for i in range(r_low):
            CNOT | (a[i], r_a[i])
        for i in range(n//2):
            CNOT | (a[r_low + i], r_a[i])
        for i in range(r_low):
            CNOT | (b[i], r_b[i])
        for i in range(n//2):
            CNOT | (b[r_low + i], r_b[i])

    # upper-part setting
    if(r_low == 1):
        c = eng.allocate_qureg(3)
        Toffoli_gate(eng, a[0], b[0], c[0])
        Toffoli_gate(eng, a[1], b[1], c[2])
        CNOT | (c[0], c[1])
        CNOT | (c[2], c[1])
        Toffoli_gate(eng, r_a, r_b, c[1])

        Uncompute(eng)
        return c, count, ancilla

    c_a = []
    c_b = []
    c_r = []

    c_a, count, ancilla = recursive_karatsuba(eng, a[0:r_low], b[0:r_low], r_low, count, ancilla)# 2 qubits     # 0~2
    c_b, count, ancilla = recursive_karatsuba(eng, a[r_low:n], b[r_low:n], n//2, count, ancilla)#2 qubits        # 3~5
    c_r, count, ancilla = recursive_karatsuba(eng, r_a[0:r_low], r_b[0:r_low], r_low, count, ancilla) #2qubits  # 6~8

    Uncompute(eng)
    #print('check initialize')
    #print_state(eng, r_a, r_low)
    result = []
    result = combine(eng, c_a, c_b, c_r, n)

    return result, count, ancilla

def Reduction(eng, result):
    CNOT | (result[8], result[0])
    CNOT | (result[12], result[0])
    CNOT | (result[13], result[0])

    CNOT | (result[12], result[1])
    CNOT | (result[13], result[2])
    CNOT | (result[8], result[3])

    CNOT | (result[11], result[4])
    CNOT | (result[12], result[5])
    CNOT | (result[11], result[7])

    #######
    CNOT | (result[14], result[8])  # 8 = green
    CNOT | (result[9], result[8])
    CNOT | (result[10], result[9])  # 9 = yellow

    CNOT | (result[12], result[14])  # 14 = blue
    CNOT | (result[13], result[10])  # 10 = orange
    CNOT | (result[11], result[10])

    #######
    CNOT | (result[14], result[3])
    CNOT | (result[14], result[7])
    CNOT | (result[10], result[3])

    CNOT | (result[10], result[6])
    CNOT | (result[8], result[1])
    CNOT | (result[8], result[4])

    CNOT | (result[9], result[2])
    CNOT | (result[9], result[5])

    return result[0:8]


def copy(eng, a, b, n):
    for i in range(n):
        CNOT | (a[i], b[i])


def Squaring(eng, vector, n):
    out = []
    CNOT | (vector[5], vector[3])
    CNOT | (vector[7], vector[2])
    CNOT | (vector[4], vector[2])
    CNOT | (vector[4], vector[0])
    CNOT | (vector[6], vector[7])
    CNOT | (vector[6], vector[0])
    CNOT | (vector[5], vector[6])
    CNOT | (vector[5], vector[1])
    CNOT | (vector[7], vector[4])
    CNOT | (vector[4], vector[5])

    out.append(vector[0])
    out.append(vector[4])
    out.append(vector[1])
    out.append(vector[5])

    out.append(vector[2])
    out.append(vector[6])
    out.append(vector[3])
    out.append(vector[7])

    return out


def combine(eng, a, b, r, n):
    if (n % 2 != 0):
        # n = 13########
        for i in range(n):
            CNOT | (a[i], r[i])
        for i in range(n - 2):
            CNOT | (b[i], r[i])

        for i in range(n // 2):
            CNOT | (a[n // 2 + 1 + i], r[i])
        for i in range(n // 2):
            CNOT | (b[i], r[n // 2 + 1 + i])

        out = []
        for i in range(n // 2 + 1):
            out.append(a[i])
        for i in range(n):
            out.append(r[i])
        for i in range((2 * n - 1) - n // 2 - 1 - n):
            out.append(b[n // 2 + i])

        return out

    half_n = int(n / 2)  # n=4
    for i in range(n - 1):
        CNOT | (a[i], r[i])
        CNOT | (b[i], r[i])
    for i in range(half_n - 1):
        CNOT | (a[half_n + i], r[i])
        CNOT | (b[i], r[half_n + i])

    result = []
    for i in range(half_n):
        result.append(a[i])
    for i in range(n - 1):
        result.append(r[i])
    for i in range(half_n):
        result.append(b[half_n - 1 + i])

    return result


def room(eng, length):
    room = eng.allocate_qureg(length)

    return room


def Toffoli_gate(eng, a, b, c):
    if (NCT):
        Toffoli | (a, b, c)
    else:
        if (resource_check):
            if (AND_check):
                ancilla = eng.allocate_qubit()
                H | c
                CNOT | (b, ancilla)
                CNOT | (c, a)
                CNOT | (c, b)
                CNOT | (a, ancilla)
                Tdag | a
                Tdag | b
                T | c
                T | ancilla
                CNOT | (a, ancilla)
                CNOT | (c, b)
                CNOT | (c, a)
                CNOT | (b, ancilla)
                H | c
                S | c

            else:
                Tdag | a
                Tdag | b
                H | c
                CNOT | (c, a)
                T | a
                CNOT | (b, c)
                CNOT | (b, a)
                T | c
                Tdag | a
                CNOT | (b, c)
                CNOT | (c, a)
                T | a
                Tdag | c
                CNOT | (b, a)
                H | c


def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def print_hex(eng, qubits, len):

    for i in reversed(range(len)):
        temp = 0
        temp = temp + int(qubits[4 * i + 3]) * 8
        temp = temp + int(qubits[4 * i + 2]) * 4
        temp = temp + int(qubits[4 * i + 1]) * 2
        temp = temp + int(qubits[4 * i])

        temp = hex(temp)
        y = temp.replace("0x", "")
        print(y, end='')

def print_state(eng, b, n):  # n = /4
    All(Measure) | b
    print('0x', end='')
    print_hex(eng, b, n)
    print('\n')

# def print_state(eng, b, n):
#     All(Measure) | b
#     print('Result : ', end='')
#     for i in range(n):
#         print(int(b[n-1-i]), end='')
#     print('\n')

global resource_check
global AND_check
global NCT

resource_check = 0
NCT =1
AND_check =0
classic = ClassicalSimulator()
eng = MainEngine(classic)
print("#### Ciphertext ####")
main(eng)
eng.flush()
print('\n')

# 분해 전
#NCT = 1
# resource_check = 1
# AND_check = 0
# Resource = ResourceCounter()
# eng = MainEngine(Resource)
# print("#### Quantum Resources costs (NCT Level) ####")
# main(eng)
# print(Resource)
# eng.flush()
# print('\n')
#
# # 분해 후
# NCT = 0
# Resource = ResourceCounter()
# eng = MainEngine(Resource)
# print("#### Decomposed Quantum Resources costs ####")
# main(eng)
# print(Resource)
# print('\n')