import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import ResourceCounter, ClassicalSimulator
from projectq.meta import Compute, Uncompute, Dagger


def Sbox_inverse(eng, u_in, q, s):
    temp = eng.allocate_qureg(8)

    X | u_in[6]
    X | u_in[5]
    X | u_in[1]
    X | u_in[0]


    # inverse L part
    u = []

    # u.append(u_in[4])
    # u.append(u_in[6])
    # u.append(u_in[5])
    # u.append(u_in[1])
    # u.append(u_in[0])
    # u.append(u_in[7])
    # u.append(u_in[3])
    # u.append(u_in[2])

    u.append(u_in[2])
    u.append(u_in[3])
    u.append(u_in[7])
    u.append(u_in[0])
    u.append(u_in[1])
    u.append(u_in[5])
    u.append(u_in[6])
    u.append(u_in[4])

    CNOT | (u[4], u[7])
    CNOT | (u[5], u[0])
    CNOT | (u[6], u[1])
    CNOT | (u[1], u[4])
    CNOT | (u[3], u[1])
    CNOT | (u[0], u[3])
    CNOT | (u[1], u[5])
    CNOT | (u[2], u[0])
    CNOT | (u[7], u[2])
    CNOT | (u[6], u[7])
    CNOT | (u[5], u[6])
    CNOT | (u[4], u[5])
    CNOT | (u[3], u[5])
    CNOT | (u[2], u[5])

#############################################

    CNOT | (u[6], q[4])
    CNOT | (u[5], q[4])
    CNOT | (u[4], q[0])
    CNOT | (u[7], q[0])  #
    CNOT | (u[1], u[3])
    CNOT | (u[4], q[2])
    CNOT | (u[2], q[1])
    CNOT | (u[3], q[3])  #
    CNOT | (q[0], q[3])  #
    CNOT | (q[4], q[6])  #
    CNOT | (u[7], q[1])
    CNOT | (q[4], q[7])
    CNOT | (q[3], q[7])
    CNOT | (u[2], u[6])
    CNOT | (u[2], u[5])
    CNOT | (u[1], u[7])
    CNOT | (u[2], q[2])
    CNOT | (u[7], q[8])
    CNOT | (q[2], q[8])  # limit
    CNOT | (u[6], q[9])
    CNOT | (q[3], q[9])
    CNOT | (u[3], u[6])
    CNOT | (u[0], q[5])
    CNOT | (q[3], q[5])
    CNOT | (u[5], u[3])
    CNOT | (u[3], q[10])
    CNOT | (u[0], u[4])
    CNOT | (q[4], u[4])
    CNOT | (q[0], q[11])
    CNOT | (u[4], q[11])
    CNOT | (u[0], u[1])
    CNOT | (u[0], q[6])  # down
    CNOT | (u[1], q[4])
    CNOT | (q[1], q[12])
    CNOT | (q[1], q[13])
    CNOT | (q[7], q[13])
    CNOT | (q[11], q[14])
    CNOT | (u[7], q[15])
    CNOT | (u[3], q[15])
    CNOT | (q[0], u[5])
    CNOT | (q[6], q[10])  #
    CNOT | (q[10], q[14])
    CNOT | (q[4], q[12])

    Toffoli_gate(eng, q[8], q[3], q[16])
    Toffoli_gate(eng, q[12], q[5], q[17])
    Toffoli_gate(eng, u[4], u[0], q[18])
    Toffoli_gate(eng, u[7], u[3], q[19])
    Toffoli_gate(eng, q[4], q[6], q[20])
    Toffoli_gate(eng, q[11], q[10], q[21])
    Toffoli_gate(eng, q[0], u[6], q[22])
    Toffoli_gate(eng, q[2], u[5], q[23])
    Toffoli_gate(eng, q[1], q[7], q[24])

    CNOT | (q[16], q[9])
    CNOT | (q[18], q[16])
    CNOT | (q[19], q[15])
    CNOT | (q[19], q[21])
    CNOT | (q[22], q[23])
    CNOT | (q[22], q[24])
    CNOT | (q[17], q[9])
    CNOT | (q[13], q[16])
    CNOT | (q[20], q[15])
    CNOT | (q[24], q[21])
    CNOT | (q[23], q[9])
    CNOT | (q[24], q[16])
    CNOT | (q[23], q[15])
    CNOT | (q[14], q[21])  # can be moved to top
    CNOT | (q[15], q[25])
    CNOT | (q[21], q[25])
    CNOT | (q[15], q[60])
    CNOT | (q[9], q[61])

    CNOT | (q[16], q[28])
    CNOT | (q[9], q[28])
    CNOT | (q[28], q[62])
    CNOT | (q[28], q[34])

    Toffoli_gate(eng, q[15], q[9], q[26])
    Toffoli_gate(eng, q[61], q[21], q[32])
    Toffoli_gate(eng, q[16], q[60], q[35])

    CNOT | (q[25], q[63])
    CNOT | (q[16], q[27])
    CNOT | (q[26], q[27])
    CNOT | (q[21], q[29])
    CNOT | (q[26], q[29])

    CNOT | (q[26], q[34])

    Toffoli_gate(eng, q[29], q[28], q[16])
    Toffoli_gate(eng, q[27], q[25], q[21])
    Toffoli_gate(eng, q[62], q[32], q[33])
    Toffoli_gate(eng, q[63], q[35], q[36])

    CNOT | (q[25], q[26])
    # CNOT | (q[30], q[16])
    CNOT | (q[34], q[33])
    # CNOT | (q[31], q[21])
    CNOT | (q[26], q[36])
    CNOT | (q[36], q[65])  #
    CNOT2(eng, q[33], q[36], q[37])
    CNOT2(eng, q[16], q[21], q[38])
    CNOT2(eng, q[16], q[33], q[39])
    CNOT2(eng, q[21], q[36], q[40])
    CNOT2(eng, q[38], q[37], q[41])
    CNOT | (q[40], q[64])

    CNOT | (q[21], q[66])
    CNOT | (q[39], q[67])
    CNOT | (q[33], q[68])
    CNOT | (q[16], q[69])
    CNOT | (q[38], q[70])
    CNOT | (q[41], q[71])
    CNOT | (q[37], q[72])

    Toffoli_gate(eng, q[40], q[3], q[42])
    Toffoli_gate(eng, q[36], q[5], q[63])
    Toffoli_gate(eng, q[21], u[0], q[44])
    Toffoli_gate(eng, q[39], u[3], q[45])
    Toffoli_gate(eng, q[33], q[6], q[46])
    Toffoli_gate(eng, q[16], q[10], q[47])
    Toffoli_gate(eng, q[38], u[6], q[48])
    Toffoli_gate(eng, q[41], u[5], q[49])
    Toffoli_gate(eng, q[37], q[7], q[50])
    Toffoli_gate(eng, q[64], q[8], q[51])
    Toffoli_gate(eng, q[65], q[12], q[52])
    # Toffoli_gate(eng, q[66], u[4], q[53], resource_check)
    Toffoli_gate(eng, q[67], u[7], q[54])
    Toffoli_gate(eng, q[68], q[4], q[55])
    Toffoli_gate(eng, q[69], q[11], q[56])
    Toffoli_gate(eng, q[70], q[0], q[57])
    Toffoli_gate(eng, q[71], q[2], q[58])
    # Toffoli_gate(eng, q[72], q[1], q[59], resource_check)


    Toffoli_gate(eng, q[66], u[4], s[2])
    Toffoli_gate(eng, q[72], q[1], s[5])

    CNOT | (q[15], q[60])
    CNOT | (q[9], q[61])
    CNOT | (q[28], q[62])
    CNOT | (q[25], q[63])
    CNOT | (q[40], q[64])
    CNOT | (q[36], q[65])
    CNOT | (q[21], q[66])
    CNOT | (q[39], q[67])
    CNOT | (q[33], q[68])
    CNOT | (q[16], q[69])
    CNOT | (q[38], q[70])
    CNOT | (q[41], q[71])
    CNOT | (q[37], q[72])
    CNOT2(eng, q[57], q[58], q[60])
    CNOT2(eng, q[46], q[52], q[61])
    CNOT2(eng, q[42], q[44], q[62])
    CNOT | (q[51], q[63])  #
    CNOT2(eng, q[50], q[54], q[64])
    CNOT | (q[45], q[65])  #
    CNOT | (q[57], q[65])  #
    CNOT | (q[58], q[66])  #
    CNOT | (q[65], q[66])  #

    # CNOT | (q[43], q[63])  # Here

    CNOT2(eng, q[42], q[63], q[67])
    CNOT2(eng, q[47], q[55], q[68])
    CNOT2(eng, q[48], q[49], q[69])
    CNOT2(eng, q[49], q[64], q[70])
    CNOT2(eng, q[56], q[62], q[71])
    CNOT2(eng, q[44], q[47], q[72])

    CNOT | (q[66], s[2])
    CNOT | (q[70], s[2])

    CNOT | (q[60], q[46])
    CNOT | (q[57], q[48])
    CNOT | (q[61], q[51])
    CNOT | (q[60], q[52])

    CNOT | (q[68], q[54])

    CNOT2(eng, q[61], q[67], s[4])

    CNOT | (q[61], q[60])

    # X | s[6]
    # X | s[5]
    # X | s[1]
    # X | s[0]

    CNOT | (q[61], s[2])
    CNOT | (q[52], s[6])

    CNOT2(eng, q[63], q[72], s[3])
    CNOT2(eng, q[54], q[62], s[0])
    CNOT2(eng, q[51], q[69], s[7])
    CNOT2(eng, q[67], q[69], s[6])
    CNOT2(eng, q[68], q[70], s[1])
    CNOT2(eng, q[71], q[48], s[5])
    # CNOT2(eng, q[71], q[53], s[2])
    CNOT | (q[71], s[2])
    CNOT | (q[64], s[5])
    CNOT | (q[66], s[7])
    # CNOT | (q[59], s[5])
    CNOT | (q[66], s[4])
    CNOT | (q[60], s[3])
    CNOT | (q[46], s[1])
    CNOT | (q[66], s[0])

    # if (round != 9):
    #     Uncompute(eng)
    #     Uncompute(eng)
    #     Uncompute(eng)
    #     Uncompute(eng)


#############################################


    s_out = []
    # s_out.append(s[4])
    # s_out.append(s[6])
    # s_out.append(s[5])
    # s_out.append(s[1])
    # s_out.append(s[0])
    # s_out.append(s[7])
    # s_out.append(s[3])
    # s_out.append(s[2])

    s_out.append(s[2])
    s_out.append(s[3])
    s_out.append(s[7])
    s_out.append(s[0])
    s_out.append(s[1])
    s_out.append(s[5])
    s_out.append(s[6])
    s_out.append(s[4])

    CNOT | (s_out[4], s_out[7])
    CNOT | (s_out[5], s_out[0])
    CNOT | (s_out[6], s_out[1])
    CNOT | (s_out[1], s_out[4])
    CNOT | (s_out[3], s_out[1])
    CNOT | (s_out[0], s_out[3])
    CNOT | (s_out[1], s_out[5])
    CNOT | (s_out[2], s_out[0])
    CNOT | (s_out[7], s_out[2])
    CNOT | (s_out[6], s_out[7])
    CNOT | (s_out[5], s_out[6])
    CNOT | (s_out[4], s_out[5])
    CNOT | (s_out[3], s_out[5])
    CNOT | (s_out[2], s_out[5])

    return s_out


def CNOT2(eng, a, b, c):
    CNOT | (a, c)
    CNOT | (b, c)


def CNOT3(eng, a, b, c, d):  # output:d
    CNOT | (a, d)
    CNOT | (b, d)
    CNOT | (c, d)


def CNOT4(eng, a, b, c, d, e):  # output:e
    CNOT3(eng, a, b, c, e)
    CNOT | (d, e)


def CNOT5(eng, a, b, c, d, e, f):  # output:f
    CNOT4(eng, a, b, c, d, f)
    CNOT | (e, f)


def CNOT6(eng, a, b, c, d, e, f, g):  # output:g
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
    q = eng.allocate_qureg(73)
    s = eng.allocate_qureg(8)
    P = 0xab

    if (resource_check != 1):
        Round_constant_XOR(eng, u, P, 8)

    u = Sbox_inverse(eng, u, q, s)

    if (resource_check != 1):
        print_state(eng, u, 2)


global resource_check
global AND_check
global NCT

resource_check = 0
AND_check = 0

# 결과값 확인
classic = ClassicalSimulator()
eng = MainEngine(classic)
print("#### Ciphertext ####")
main(eng)
eng.flush()
print('\n')

resource_check = 1
NCT = 0
AND_check = 0
Resource = ResourceCounter()
eng = MainEngine(Resource)
main(eng)
print('\n')
print(Resource)