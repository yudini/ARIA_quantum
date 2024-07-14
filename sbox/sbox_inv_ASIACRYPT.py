import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import ResourceCounter, ClassicalSimulator
from projectq.meta import Compute, Uncompute, Dagger

def Sbox_inverse(eng, input, t, s):

    temp = eng.allocate_qureg(8)

    u_in =[]

    for i in range(8):
        u_in.append(input[7-i])


    X | u_in[1]
    X | u_in[2]
    X | u_in[6]
    X | u_in[7]

    # inverse L part
    u=[]
    u.append(u_in[3])
    u.append(u_in[1])
    u.append(u_in[2])
    u.append(u_in[6])
    u.append(u_in[7])
    u.append(u_in[0])
    u.append(u_in[4])
    u.append(u_in[5])

    CNOT | (u[3], u[0])
    CNOT | (u[2], u[7])
    CNOT | (u[1], u[6])
    CNOT | (u[6], u[3])
    CNOT | (u[4], u[6])
    CNOT | (u[7], u[4])
    CNOT | (u[6], u[2])
    CNOT | (u[5], u[7])
    CNOT | (u[0], u[5])
    CNOT | (u[1], u[0])
    CNOT | (u[2], u[1])
    CNOT | (u[3], u[2])
    CNOT | (u[4], u[2])
    CNOT | (u[5], u[2])


    #####################

    CNOT | (u[6], u[0])
    CNOT3(eng, u[2], u[4], u[5], u[6])
    Toffoli_gate(eng, u[0], u[6], t[0])
    CNOT | (t[0], t[1])

    CNOT2(eng, u[2], u[7], u[1])
    CNOT4(eng, u[1], u[4], u[5], u[6], u[2])
    Toffoli_gate(eng, u[1], u[2], t[1])

    CNOT3(eng, u[1], u[2], u[3], u[0])
    CNOT2(eng, u[1], u[7], u[6])
    Toffoli_gate(eng, u[0], u[6], t[2])
    CNOT | (t[2], t[3])

    CNOT | (u[3], u[5])
    CNOT4(eng, u[2], u[4], u[6], u[7], u[0])
    Toffoli_gate(eng, u[5], u[0], t[3])
    CNOT | (t[3], t[1])

    CNOT2(eng, u[3], u[4], u[0])
    CNOT | (u[0], t[1])

    CNOT4(eng, u[1], u[2], u[6], u[7], u[0])
    CNOT | (u[7], u[6])
    Toffoli_gate(eng, u[0], u[6], t[0])

    CNOT2(eng, u[2], u[5], u[0])
    CNOT3(eng, u[0], u[3], u[4], u[5])
    Toffoli_gate(eng, u[0], u[5], t[4])

    CNOT | (t[4], t[3])
    CNOT2(eng, u[1], u[3], u[0])
    CNOT | (u[7], u[5])
    Toffoli_gate(eng, u[0], u[5], t[3])

    CNOT4(eng, u[2], u[4], u[5], u[6], u[1])
    CNOT | (u[1], t[3])

    CNOT3(eng, u[2], u[4], u[6], u[1])
    CNOT | (u[2], u[0])
    Toffoli_gate(eng, u[1], u[0], t[2])
    CNOT | (t[2], t[0])
    CNOT | (t[4], t[2])

    CNOT4(eng, u[1], u[2], u[3], u[5], u[0])
    CNOT | (u[7], u[5])
    Toffoli_gate(eng, u[0], u[5], t[4])

    CNOT4(eng, u[3], u[4], u[5], u[6], u[2])
    CNOT | (u[2], t[0])

    CNOT3(eng, u[3], u[5], u[7], u[1])
    Toffoli_gate(eng, u[1], u[7], t[2])

    CNOT4(eng, u[2], u[4], u[6], u[7], u[0])
    CNOT | (u[0], t[2])

    Toffoli_gate(eng, t[1], t[3], t[4])
    CNOT | (t[0], t[1])
    CNOT | (t[2], t[4])
    Toffoli_gate(eng, t[1], t[4], t[5])
    CNOT | (t[0], t[5])
    CNOT2(eng, t[0], t[5], t[1])
    CNOT2(eng, t[2], t[0], t[4])
    CNOT | (t[4], t[5])
    Toffoli_gate(eng, t[0], t[5], t[1])
    CNOT | (t[2], t[3])
    Toffoli_gate(eng, t[3], t[4], t[2])
    Toffoli_gate(eng, t[0], t[5], t[4])
    Toffoli_gate(eng, t[2], t[4], t[3])
    Toffoli_gate(eng, t[0], t[5], t[4])
    CNOT | (t[4], t[5])

    CNOT4(eng, u[2], u[4], u[5], u[6], u[3])
    CNOT3(eng, u[0], u[3], u[7], u[4])
    CNOT | (u[6], u[2])
    CNOT | (u[7], u[5])
    CNOT5(eng, u[0], u[1], u[3], u[4], u[5], u[6])

    ###############
    # R2

    CNOT3(eng, u[0], u[4], u[2], u[1])
    CNOT | (t[1], t[3])

    Toffoli_gate(eng, t[3], u[1], s[5])
    CNOT | (s[5], s[6])

    CNOT3(eng, u[0], u[4], u[2], u[1])
    CNOT | (t[1], t[3])

    Toffoli_gate(eng, t[2], u[2], s[2])
    CNOT | (s[2], s[5])
    CNOT | (s[6], s[2])

    Toffoli_gate(eng, t[1], u[5], s[4])
    CNOT | (s[4], s[1])
    CNOT | (s[4], s[3])

    CNOT | (u[7], u[5])
    CNOT | (t[5], t[1])

    Toffoli_gate(eng, t[1], u[5], s[7])
    CNOT | (s[7], s[1])
    CNOT | (s[7], s[3])
    CNOT | (s[7], s[4])

    CNOT | (u[7], u[5])
    CNOT | (t[5], t[1])

    Toffoli_gate(eng, t[5], u[7], s[7])
    CNOT | (s[7], s[2])
    CNOT | (s[7], s[5])
    CNOT | (s[7], s[6])

    CNOT5(eng, u[3], u[0], u[4], u[5], u[6], u[1])

    Toffoli_gate(eng, t[2], u[1], s[7])
    CNOT | (s[7], s[2])
    CNOT | (s[7], s[4])
    CNOT | (s[7], s[5])

    CNOT5(eng, u[3], u[0], u[4], u[5], u[6], u[1])

    CNOT | (u[3], u[2])
    CNOT | (t[2], t[3])

    Toffoli_gate(eng, t[3], u[2], s[7])
    CNOT | (s[7], s[2])
    CNOT | (s[7], s[5])

    CNOT | (u[3], u[2])
    CNOT | (t[2], t[3])

    Toffoli_gate(eng, t[3], u[3], s[7])
    CNOT | (s[7], s[6])

    CNOT2(eng, u[3], u[2], u[6])
    CNOT | (t[3], t[2])

    Toffoli_gate(eng, t[2], u[6], s[0])
    CNOT | (s[0], s[4])
    CNOT | (s[0], s[6])
    CNOT | (s[0], s[7])

    CNOT2(eng, u[3], u[2], u[6])
    CNOT | (t[3], t[2])

    CNOT4(eng, u[0], u[4], u[2], u[5], u[1])

    Toffoli_gate(eng, t[3], u[1], s[0])
    CNOT | (s[0], s[1])
    CNOT | (s[0], s[2])
    CNOT | (s[0], s[3])
    CNOT | (s[0], s[4])
    CNOT | (s[0], s[5])
    CNOT | (s[0], s[6])

    CNOT4(eng, u[0], u[4], u[2], u[5], u[1])

    CNOT6(eng, u[7], u[0], u[4], u[5], u[6], u[1], u[3])
    CNOT | (t[2], t[5])

    Toffoli_gate(eng, t[5], u[3], s[0])
    CNOT | (s[0], s[2])
    CNOT | (s[0], s[5])
    CNOT | (s[0], s[6])

    CNOT6(eng, u[7], u[0], u[4], u[5], u[6], u[1], u[3])
    CNOT | (t[2], t[5])

    CNOT4(eng, u[7], u[2], u[5], u[6], u[3])
    CNOT3(eng, t[2], t[1], t[3], t[5])

    Toffoli_gate(eng, t[5], u[3], s[0])
    CNOT | (s[0], s[3])
    CNOT | (s[0], s[4])
    CNOT | (s[0], s[5])
    CNOT | (s[0], s[6])

    CNOT4(eng, u[7], u[2], u[5], u[6], u[3])
    CNOT3(eng, t[2], t[1], t[3], t[5])

    CNOT2(eng, u[3], u[4], u[2])
    CNOT | (t[1], t[5])

    Toffoli_gate(eng, t[5], u[2], s[0])
    CNOT | (s[0], s[5])

    CNOT2(eng, u[3], u[4], u[2])
    CNOT | (t[1], t[5])

    CNOT3(eng, u[3], u[4], u[2], u[1])

    Toffoli_gate(eng, t[1], u[1], s[0])
    CNOT | (s[0], s[2])
    CNOT | (s[0], s[6])
    CNOT | (s[0], s[7])

    CNOT3(eng, u[3], u[4], u[2], u[1])

    CNOT | (u[2], u[1])
    CNOT | (t[2], t[5])

    Toffoli_gate(eng, t[5], u[1], s[0])
    CNOT | (s[0], s[2])

    CNOT | (u[2], u[1])
    CNOT | (t[2], t[5])

    CNOT3(eng, t[2], t[1], t[3], t[5])

    Toffoli_gate(eng, t[5], u[4], s[0])
    CNOT | (s[0], s[1])
    CNOT | (s[0], s[3])
    CNOT | (s[0], s[4])
    CNOT | (s[0], s[5])
    CNOT | (s[0], s[6])
    CNOT | (s[0], s[7])

    CNOT3(eng, t[2], t[1], t[3], t[5])

    CNOT2(eng, u[4], u[2], u[1])
    CNOT | (t[1], t[3])

    Toffoli_gate(eng, t[3], u[1], s[2])

    CNOT2(eng, u[4], u[2], u[1])
    CNOT | (t[1], t[3])


    Toffoli_gate(eng, t[5], u[1], s[5])


    s_out = []
    s_out.append(s[3])
    s_out.append(s[1])
    s_out.append(s[2])
    s_out.append(s[6])
    s_out.append(s[7])
    s_out.append(s[0])
    s_out.append(s[4])
    s_out.append(s[5])

    CNOT | (s_out[3], s_out[0])
    CNOT | (s_out[2], s_out[7])
    CNOT | (s_out[1], s_out[6])
    CNOT | (s_out[6], s_out[3])
    CNOT | (s_out[4], s_out[6])
    CNOT | (s_out[7], s_out[4])
    CNOT | (s_out[6], s_out[2])
    CNOT | (s_out[5], s_out[7])
    CNOT | (s_out[0], s_out[5])
    CNOT | (s_out[1], s_out[0])
    CNOT | (s_out[2], s_out[1])
    CNOT | (s_out[3], s_out[2])
    CNOT | (s_out[4], s_out[2])
    CNOT | (s_out[5], s_out[2])

    output = []

    for i in range(8):
        output.append(s_out[7-i])

    return output

def CNOT2(eng, a, b, c):
    CNOT | (a, c)
    CNOT | (b, c)

def CNOT3(eng, a, b, c, d): #output:d
    CNOT |(a,d)
    CNOT | (b,d)
    CNOT | (c,d)

def CNOT4(eng, a, b, c, d, e): #output:e
    CNOT3(eng, a, b, c, e)
    CNOT | (d,e)

def CNOT5(eng, a, b, c, d, e, f): #output:f
    CNOT4(eng, a, b, c, d, f)
    CNOT | (e, f)

def CNOT6(eng, a, b, c, d, e, f, g): #output:g
    CNOT5(eng, a, b, c, d, e, g)
    CNOT | (f, g)


def Toffoli_gate(eng, a, b, c):

    Toffoli | (a,b,c)

    # if (resource_check):
    #     Tdag | a
    #     Tdag | b
    #     H | c
    #     CNOT | (c, a)
    #     T | a
    #     CNOT | (b, c)
    #     CNOT | (b, a)
    #     T | c
    #     Tdag | a
    #     CNOT | (b, c)
    #     CNOT | (c, a)
    #     T | a
    #     Tdag | c
    #     CNOT | (b, a)
    #     H | c
    # else:
    #     Toffoli | (a, b, c)

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

def main(eng):
    u = eng.allocate_qureg(8)
    t = eng.allocate_qureg(6)
    s = eng.allocate_qureg(8)
    P = 0x5e

    if (resource_check != 1):
        Round_constant_XOR(eng, u, P, 8)

    u = Sbox_inverse(eng,u,t,s)

    if(resource_check!=1):
        print_state(eng,u,2)

global resource_check
global AND_check
global NCT

resource_check = 0
AND_check = 0

# 결과값 확인
# classic = ClassicalSimulator()
# eng = MainEngine(classic)
# print("#### Ciphertext ####")
# main(eng)
# eng.flush()
# print('\n')


resource_check = 1
NCT = 1
AND_check = 0
Resource = ResourceCounter()
eng = MainEngine(Resource)
main(eng)
print('\n')
print(Resource)