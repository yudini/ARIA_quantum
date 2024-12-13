
import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control, Dagger


def ARIA256(eng):
    r = 16  # round number
    bit = 8
    n = 128

    CK_1 = 0x517cc1b727220a94fe13abe8fa9a6ee0
    CK_2 = 0x6db14acc9e21c820ff28b1d5ef5de2b0
    CK_3 = 0xdb92371d2126e9700324977504e8c90e

    PT = eng.allocate_qureg(n)
    w0 = eng.allocate_qureg(n)
    RK = eng.allocate_qureg(n)

    w1 = eng.allocate_qureg(n)
    w2 = eng.allocate_qureg(n)
    w3 = eng.allocate_qureg(n)

    #ancilla = eng.allocate_qureg(304)  # 38*4*4 (병렬)

    vect_a = 0x63
    vect_b = 0xe2

    div = pow(16, 32)
    KL = MK // div
    KR = MK % div

    if (resource_check != 1):
        Round_constant_XOR(eng, w0, KL, n)
        Round_constant_XOR(eng, PT, P, n)
        Round_constant_XOR(eng, w1, CK_3, n)
        Round_constant_XOR(eng, w2, CK_1, n)
        Round_constant_XOR(eng, w3, CK_2, n)

    w1 = RoundOdd(eng, w1, w0, vect_a, vect_b)

    if (resource_check != 1):
        print_state(eng, w1, 32)
        print_state(eng, w2, 32)
        print_state(eng, w3, 32)

    Round_constant_XOR(eng, w1, KR, n)  # 256-128 = 128

    ## W2 계산
    w2 = RoundEven(eng, w2, w1, vect_a, vect_b)

    for i in range(n):
        CNOT | (w0[i], w2[i])

    if (resource_check != 1):
        print_state(eng, w1, 32)
        print_state(eng, w2, 32)
        print_state(eng, w3, 32)

    ## w3 계산
    w3 = RoundOdd(eng, w3, w2, vect_a, vect_b)


    if (resource_check != 1):
        print_state(eng, w1, 32)
        print_state(eng, w2, 32)
        print_state(eng, w3, 32)

    for i in range(n):
        CNOT | (w1[i], w3[i])

    W = [w0, w1, w2, w3]
    Wnum = [19, 31, 67, 97, 109]

    #### Round Operation ###
    for i in range(r):
        print("Round {0}".format(i + 1))

        #### Generate round-key ####
        key_generation(eng, KL, W[(i + 1) % 4], RK, Wnum[i // 4], W[i % 4], i, n)

        ## Final round ##
        if i == r - 1:
            PT = RoundFinal(eng, PT, RK, vect_a, vect_b)

        ## odd round ##
        elif i % 2 + 1 == 1:
            PT = RoundOdd(eng, PT, RK, vect_a, vect_b)

        ## even round ##
        elif i % 2 == 1:
            PT = RoundEven(eng, PT, RK, vect_a, vect_b)

        rev_key_generation(eng, KL, W[(i + 1) % 4], RK, Wnum[i // 4], W[i % 4], i, n)

    #### Final Round ####
    key_generation(eng, KL, W[(r + 1) % 4], RK, Wnum[r // 4], W[r % 4], r, n)

    for i in range(n):
        CNOT | (RK[i], PT[i])

    #### print Ciphertext ####
    if (resource_check != 1):
        print("### Ciphertext ###")
        print_state(eng, PT, n // 4) # 58a875e6044ad7fffa4f58420f7f442d


def key_generation(eng, KL, w_rot, rk, shift, w_i, cnt, n):
    ## Initialization ##
    if cnt % 4 == 0:
        Round_constant_XOR(eng, rk, KL, n)  # RK1,5,9,13,17 = w0
    else:
        for i in range(128):
            CNOT | (w_i[i], rk[i])

    ## Generation ##
    for i in range(shift, n):
        CNOT | (w_rot[i], rk[i - shift])
    for i in range(0, shift):
        CNOT | (w_rot[i], rk[n - shift + i])


def rev_key_generation(eng, KL, w_rot, rk, shift, w_i, cnt, n):
    with Dagger(eng):
        ## Initialization ##
        if cnt % 4 == 0:
            Round_constant_XOR(eng, rk, KL, n)  # RK1,5,9,13 = w0
        else:
            for i in range(128):
                CNOT | (w_i[i], rk[i])

        ## Generation ##
        for i in range(shift, n):
            CNOT | (w_rot[i], rk[i - shift])
        for i in range(0, shift):
            CNOT | (w_rot[i], rk[n - shift + i])


def RoundOdd(eng, input, key, vect_a, vect_b):
    ## input ^ key
    for i in range(128):
        CNOT | (key[i], input[i])

    ## S-box (Type1)
    input[120:] = Sbox(eng, input[120:])
    input[112:120] = Sbox2(eng, input[112:120], vect_a, vect_b)
    input[104:112] = SboxInv(eng, input[104:112], vect_a)
    input[96:104] = Sbox2Inv(eng, input[96:104], vect_a, vect_b)

    input[88:96] = Sbox(eng, input[88:96])
    input[80:88] = Sbox2(eng, input[80:88], vect_a, vect_b)
    input[72:80] = SboxInv(eng, input[72:80], vect_a)
    input[64:72] = Sbox2Inv(eng, input[64:72], vect_a, vect_b)

    input[56:64] = Sbox(eng, input[56:64])
    input[48:56] = Sbox2(eng, input[48:56], vect_a, vect_b)
    input[40:48] = SboxInv(eng, input[40:48], vect_a)
    input[32:40] = Sbox2Inv(eng, input[32:40], vect_a, vect_b)

    input[24:32] = Sbox(eng, input[24:32])
    input[16:24] = Sbox2(eng, input[16:24], vect_a, vect_b)
    input[8:16] = SboxInv(eng, input[8:16], vect_a)
    input[0:8] = Sbox2Inv(eng, input[0:8], vect_a, vect_b)
    # input[120:] = Sbox1(eng, input[120:])
    # input[112:120] = Sbox2(eng, input[112:120], vect_b, bit, ancilla[:38])
    # input[104:112] = Sbox1Inv(eng, input[104:112])
    # input[96:104] = Sbox2Inv(eng, input[96:104], vect_b, bit, ancilla[38:76])
    #
    # input[88:96] = Sbox1(eng, input[88:96])
    # input[80:88] = Sbox2(eng, input[80:88], vect_b, bit, ancilla[76:114])
    # input[72:80] = Sbox1Inv(eng, input[72:80])
    # input[64:72] = Sbox2Inv(eng, input[64:72], vect_b, bit, ancilla[114:152])
    #
    # input[56:64] = Sbox1(eng, input[56:64])
    # input[48:56] = Sbox2(eng, input[48:56], vect_b, bit, ancilla[152:190])
    # input[40:48] = Sbox1Inv(eng, input[40:48])
    # input[32:40] = Sbox2Inv(eng, input[32:40], vect_b, bit, ancilla[190:228])
    #
    # input[24:32] = Sbox1(eng, input[24:32])
    # input[16:24] = Sbox2(eng, input[16:24], vect_b, bit, ancilla[228:266])
    # input[8:16] = Sbox1Inv(eng, input[8:16])
    # input[0:8] = Sbox2Inv(eng, input[0:8], vect_b, bit, ancilla[266:])

    ## Diffusion layer
    input = Diffusion(eng, input)

    return input


def RoundEven(eng, input, key, vect_a, vect_b):
    ## input ^ key
    for i in range(128):
        CNOT | (key[i], input[i])
    # eng, u, q, s, flag, round, resource_check
    ## S-box type2
    input[120:] = SboxInv(eng, input[120:], vect_a)
    input[112:120] = Sbox2Inv(eng, input[112:120], vect_a, vect_b)
    input[104:112] = Sbox(eng, input[104:112])
    input[96:104] = Sbox2(eng, input[96:104], vect_a, vect_b)

    input[88:96] = SboxInv(eng, input[88:96], vect_a)
    input[80:88] = Sbox2Inv(eng, input[80:88], vect_a, vect_b)
    input[72:80] = Sbox(eng, input[72:80])
    input[64:72] = Sbox2(eng, input[64:72], vect_a, vect_b)

    input[56:64] = SboxInv(eng, input[56:64], vect_a)
    input[48:56] = Sbox2Inv(eng, input[48:56], vect_a, vect_b)
    input[40:48] = Sbox(eng, input[40:48])
    input[32:40] = Sbox2(eng, input[32:40], vect_a, vect_b)

    input[24:32] = SboxInv(eng, input[24:32], vect_a)
    input[16:24] = Sbox2Inv(eng, input[16:24], vect_a, vect_b)
    input[8:16] = Sbox(eng, input[8:16])
    input[0:8] = Sbox2(eng, input[0:8], vect_a, vect_b)
    # input[120:] = Sbox1Inv(eng, input[120:])
    # input[112:120] = Sbox2Inv(eng, input[112:120], vect_b, bit, ancilla[:38])
    # input[104:112] = Sbox1(eng, input[104:112])
    # input[96:104] = Sbox2(eng, input[96:104], vect_b, bit, ancilla[38:76])
    #
    # input[88:96] = Sbox1Inv(eng, input[88:96])
    # input[80:88] = Sbox2Inv(eng, input[80:88], vect_b, bit, ancilla[76:114])
    # input[72:80] = Sbox1(eng, input[72:80])
    # input[64:72] = Sbox2(eng, input[64:72], vect_b, bit, ancilla[114:152])
    #
    # input[56:64] = Sbox1Inv(eng, input[56:64])
    # input[48:56] = Sbox2Inv(eng, input[48:56], vect_b, bit, ancilla[152:190])
    # input[40:48] = Sbox1(eng, input[40:48])
    # input[32:40] = Sbox2(eng, input[32:40], vect_b, bit, ancilla[190:228])
    #
    # input[24:32] = Sbox1Inv(eng, input[24:32])
    # input[16:24] = Sbox2Inv(eng, input[16:24], vect_b, bit, ancilla[228:266])
    # input[8:16] = Sbox1(eng, input[8:16])
    # input[0:8] = Sbox2(eng, input[0:8], vect_b, bit, ancilla[266:])

    ## Diffusion layer
    input = Diffusion(eng, input)

    return input


def RoundFinal(eng, input, key, vect_a, vect_b):
    ## input ^ key
    for i in range(128):
        CNOT | (key[i], input[i])

    ## S-box type2
    input[120:] = SboxInv(eng, input[120:], vect_a)
    input[112:120] = Sbox2Inv(eng, input[112:120], vect_a, vect_b)
    input[104:112] = Sbox(eng, input[104:112])
    input[96:104] = Sbox2(eng, input[96:104], vect_a, vect_b)

    input[88:96] = SboxInv(eng, input[88:96], vect_a)
    input[80:88] = Sbox2Inv(eng, input[80:88], vect_a, vect_b)
    input[72:80] = Sbox(eng, input[72:80])
    input[64:72] = Sbox2(eng, input[64:72], vect_a, vect_b)

    input[56:64] = SboxInv(eng, input[56:64], vect_a)
    input[48:56] = Sbox2Inv(eng, input[48:56], vect_a, vect_b)
    input[40:48] = Sbox(eng, input[40:48])
    input[32:40] = Sbox2(eng, input[32:40], vect_a, vect_b)

    input[24:32] = SboxInv(eng, input[24:32], vect_a)
    input[16:24] = Sbox2Inv(eng, input[16:24], vect_a, vect_b)
    input[8:16] = Sbox(eng, input[8:16])
    input[0:8] = Sbox2(eng, input[0:8], vect_a, vect_b)

    return input
    # input[120:] = Sbox1Inv(eng, input[120:])
    # input[112:120] = Sbox2Inv(eng, input[112:120], vect_b, bit, ancilla[:38])
    # input[104:112] = Sbox1(eng, input[104:112])
    # input[96:104] = Sbox2(eng, input[96:104], vect_b, bit, ancilla[38:76])
    #
    # input[88:96] = Sbox1Inv(eng, input[88:96])
    # input[80:88] = Sbox2Inv(eng, input[80:88], vect_b, bit, ancilla[76:114])
    # input[72:80] = Sbox1(eng, input[72:80])
    # input[64:72] = Sbox2(eng, input[64:72], vect_b, bit, ancilla[114:152])
    #
    # input[56:64] = Sbox1Inv(eng, input[56:64])
    # input[48:56] = Sbox2Inv(eng, input[48:56], vect_b, bit, ancilla[152:190])
    # input[40:48] = Sbox1(eng, input[40:48])
    # input[32:40] = Sbox2(eng, input[32:40], vect_b, bit, ancilla[190:228])
    #
    # input[24:32] = Sbox1Inv(eng, input[24:32])
    # input[16:24] = Sbox2Inv(eng, input[16:24], vect_b, bit, ancilla[228:266])
    # input[8:16] = Sbox1(eng, input[8:16])
    # input[0:8] = Sbox2(eng, input[0:8], vect_b, bit, ancilla[266:])

    return input


def Diffusion(eng, a):
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


def CNOT8bit(a, an, b, bn):
    for i in range(8):
        CNOT | (a[127 - an * 8 - i], b[127 - bn * 8 - i])


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


def Sbox(eng, u):

    q = eng.allocate_qureg(73)
    s = eng.allocate_qureg(8)
#   with Compute(eng) :
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
    # Toffoli_gate(eng, q[66], u[4], q[53])
    Toffoli_gate(eng, q[67], u[7], q[54])
    Toffoli_gate(eng, q[68], q[4], q[55])
    Toffoli_gate(eng, q[69], q[11], q[56])
    Toffoli_gate(eng, q[70], q[0], q[57])
    Toffoli_gate(eng, q[71], q[2], q[58])
    # Toffoli_gate(eng, q[72], q[1], q[59], resource_check)


    Toffoli_gate(eng, q[66], u[4], s[2])
    Toffoli_gate(eng, q[72], q[1], s[5])

   # with Compute(eng):

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

#    with Compute(eng):

    CNOT | (q[60], q[46])
    CNOT | (q[57], q[48])
    CNOT | (q[61], q[51])
    CNOT | (q[60], q[52])

    CNOT | (q[68], q[54])

    CNOT2(eng, q[61], q[67], s[4])

    # with Compute(eng):
    #     CNOT | (q[61], q[60])

    CNOT | (q[61], q[60])

    X | s[6]
    X | s[5]
    X | s[1]
    X | s[0]

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

    return s


def SboxInv(eng, u,vect_a):
    u = MatrixProductS1Inv(eng, u, vect_a, 8)

    u = Sbox(eng, u)
    # u = reverse_L(eng, u)
    u = MatrixProductS1Inv(eng, u, vect_a, 8)
    return u


def Sbox2(eng, u,vect_a, vect_b):
    # x^-1, get inverse
    u = Sbox(eng, u)
    # u = reverse_L(eng, u)
    u = MatrixProductS1Inv(eng, u, vect_a, 8)

    # ------------------------------------
    u = MatrixProductS2(eng, u, vect_b, 8)

    return u


def Sbox2Inv(eng, u, vect_a, vect_b):
    u = MatrixProductS2Inv(eng, u, vect_b, 8)

    u = Sbox(eng, u)
    # u = reverse_L(eng, u)
    u = MatrixProductS1Inv(eng, u, vect_a, 8)

    return u



def getInverse(eng, a, a1, ancilla, n):  # n=8
    count = 0

    copy(eng, a, a1, n)  # a1에 a 복사해서 넣기
    a1 = Squaring(eng, a1, n)  # a1 = a^2

    # a2 = a * a^2
    a2 = []
    a2, count, ancilla = recursive_karatsuba(eng, a, a1, n, count, ancilla)
    a2 = Reduction(eng, a2)

    # a = a^64 -> a1 = a1 ^ 32
    a1 = Squaring(eng, a1, n)
    a1 = Squaring(eng, a1, n)
    a1 = Squaring(eng, a1, n)
    a1 = Squaring(eng, a1, n)
    a1 = Squaring(eng, a1, n)  # ^64

    # a3 = (a*a^2)*a^64
    count = 0
    a3 = []
    a3, count, ancilla = recursive_karatsuba(eng, a2, a1, n, count, ancilla)
    a3 = Reduction(eng, a3)

    # a2 = (a * a^2)^4
    a2 = Squaring(eng, a2, n)
    a2 = Squaring(eng, a2, n)

    # a4 = (a*a^2)*(a * a^2)^4*a^64
    count = 0
    a4 = []
    a4, count, ancilla = recursive_karatsuba(eng, a2, a3, n, count, ancilla)
    a4 = Reduction(eng, a4)

    # a2 = (a * a^2)^16
    a2 = Squaring(eng, a2, n)
    a2 = Squaring(eng, a2, n)

    # a5 = (a*a^2)*(a * a^2)^4*(a * a^2)^16*a^64
    count = 0
    a5 = []
    a5, count, ancilla = recursive_karatsuba(eng, a2, a4, n, count, ancilla)
    a5 = Reduction(eng, a5)

    # a5 = ((a*a^2)*(a * a^2)^4*(a * a^2)^16*a^64)^2 // a^254
    a5 = Squaring(eng, a5, n)

    return a5


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



def recursive_karatsuba(eng, a, b, n, count, ancilla):

    if (n == 1):
        c = eng.allocate_qubit()
        Toffoli_gate(eng, a, b, c)

        return c, count, ancilla

    c_len = 3 ** math.log(n, 2)
    r_low = n // 2

    if (n % 2 != 0):
        r_low = r_low + 1

    r_a = []
    r_b = []

    # Provide rooms and prepare operands
    r_a = ancilla[count:count + r_low]

    count = count + r_low

    r_b = ancilla[count:count + r_low]

    count = count + r_low

    with Compute(eng):  # r_a, r_b
        for i in range(r_low):
            CNOT | (a[i], r_a[i])
        for i in range(n // 2):
            CNOT | (a[r_low + i], r_a[i])
        for i in range(r_low):
            CNOT | (b[i], r_b[i])
        for i in range(n // 2):
            CNOT | (b[r_low + i], r_b[i])

    # upper-part setting
    if (r_low == 1):
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

    c_a, count, ancilla = recursive_karatsuba(eng, a[0:r_low], b[0:r_low], r_low, count, ancilla)  # 2 qubits     # 0~2
    c_b, count, ancilla = recursive_karatsuba(eng, a[r_low:n], b[r_low:n], n // 2, count, ancilla)  # 2 qubits    # 3~5
    c_r, count, ancilla = recursive_karatsuba(eng, r_a[0:r_low], r_b[0:r_low], r_low, count, ancilla)  # 2qubits  # 6~8

    Uncompute(eng)
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


def copy(eng, a, b, length):
    for i in range(length):
        CNOT | (a[i], b[i])


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


def print_state_a(eng, b, n):
    All(Measure) | b
    for i in range(n):
        print(int(b[n - 1 - i]), end='')
    print('\n')


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


global resource_check
global AND_check
global NCT
global MK
global P

MK = 0x00112233445566778899aabbccddeeff00112233445566778899aabbccddeeff
P = 0x11111111aaaaaaaa11111111bbbbbbbb
NCT = 1
resource_check = 0
classic = ClassicalSimulator()
eng = MainEngine(classic)
ARIA256(eng)
eng.flush()
print('\n')

# 분해 전
NCT = 1
resource_check = 1
AND_check = 0
Resource = ResourceCounter()
eng = MainEngine(Resource)
ARIA256(eng)
print('\n')
print(Resource)
eng.flush()

# 분해 후
resource_check = 1
NCT = 0
AND_check = 0
Resource = ResourceCounter()
eng = MainEngine(Resource)
ARIA256(eng)
print('\n')
print(Resource)