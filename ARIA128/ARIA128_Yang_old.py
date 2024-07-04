import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import ResourceCounter, ClassicalSimulator
from projectq.meta import Compute, Uncompute, Dagger


def ARIA128(eng):
    r = 12  # round number
    bit = 8
    n = 128

    CK_1 = 0x517cc1b727220a94fe13abe8fa9a6ee0
    CK_2 = 0x6db14acc9e21c820ff28b1d5ef5de2b0
    CK_3 = 0xdb92371d2126e9700324977504e8c90e

    PT = eng.allocate_qureg(n)
    K = eng.allocate_qureg(n)
    RK = eng.allocate_qureg(n)  # Round key
    w1 = eng.allocate_qureg(n)
    w2 = eng.allocate_qureg(n)
    w3 = eng.allocate_qureg(n)

    ancilla = eng.allocate_qureg(608)  # 38*4*4 (병렬)

    vect_a = 0x63
    vect_b = 0xe2

    if (resource_check != 1):
        Round_constant_XOR(eng, K, MK, n)
        Round_constant_XOR(eng, PT, P, n)
        Round_constant_XOR(eng, w1, CK_1, n)
        Round_constant_XOR(eng, w2, CK_2, n)
        Round_constant_XOR(eng, w3, CK_3, n)

    w0 = K[0:n]
    w1 = RoundOdd(eng, w1, w0, vect_a, vect_b, bit, ancilla)

    ## W2 계산
    w2 = RoundEven(eng, w2, w1, vect_a, vect_b, bit, ancilla)

    for i in range(n):
        CNOT | (w0[i], w2[i])

    ## w3 계산
    w3 = RoundOdd(eng, w3, w2, vect_a, vect_b, bit, ancilla)
    for i in range(n):
        CNOT | (w1[i], w3[i])

    W = [w0, w1, w2, w3]
    Wnum = [19, 31, 67, 97, 109]

    #### Round Operation ###
    for i in range(r):
        print("Round {0}".format(i + 1))

        #### Generate round-key ####
        key_generation(W[(i + 1) % 4], RK, Wnum[i // 4], W[i % 4], i, n)

        ## Final round ##
        if i == r - 1:
            PT = RoundFinal(eng, PT, RK, vect_a, vect_b, bit, ancilla)

        ## odd round ##
        elif i % 2 + 1 == 1:
            PT = RoundOdd(eng, PT, RK, vect_a, vect_b, bit, ancilla)

        ## even round ##
        elif i % 2 == 1:
            PT = RoundEven(eng, PT, RK, vect_a, vect_b, bit, ancilla)

        rev_key_generation(W[(i + 1) % 4], RK, Wnum[i // 4], W[i % 4], i, n)

    #### Final Round ####
    key_generation(W[(r + 1) % 4], RK, Wnum[r // 4], W[r % 4], r, n)
    for i in range(n):
        CNOT | (RK[i], PT[i])

    #### print Ciphertext ####
    if (resource_check != 1):
        print_state(eng, PT, n // 4)  # c6ecd08e22c30abdb215cf74e2075e6e


def key_generation(w_rot, rk, shift, w_i, cnt, n):
    ## Initialization ##
    if cnt % 4 == 0:
        Round_constant_XOR(eng, rk, MK, n)  # RK1,5,9,13 = w0
    else:
        for i in range(128):
            CNOT | (w_i[i], rk[i])

    ## Generation ##
    for i in range(shift, n):
        CNOT | (w_rot[i], rk[i - shift])
    for i in range(0, shift):
        CNOT | (w_rot[i], rk[n - shift + i])


def rev_key_generation(w_rot, rk, shift, w_i, cnt, n):
    with Dagger(eng):
        ## Initialization ##
        if cnt % 4 == 0:
            Round_constant_XOR(eng, rk, MK, n)  # RK1,5,9,13 = w0
        else:
            for i in range(128):
                CNOT | (w_i[i], rk[i])

        ## Generation ##
        for i in range(shift, n):
            CNOT | (w_rot[i], rk[i - shift])
        for i in range(0, shift):
            CNOT | (w_rot[i], rk[n - shift + i])


def RoundOdd(eng, input, key, vect_a, vect_b, bit, ancilla):
    ## input ^ key
    for i in range(128):
        CNOT | (key[i], input[i])

    ## S-box (Type1)
    input[120:] = Sbox1(eng, input[120:], vect_a, bit, ancilla[:38])
    input[112:120] = Sbox2(eng, input[112:120], vect_b, bit, ancilla[38:76])
    input[104:112] = Sbox1Inv(eng, input[104:112], vect_a, bit, ancilla[76:114])
    input[96:104] = Sbox2Inv(eng, input[96:104], vect_b, bit, ancilla[114:152])

    input[88:96] = Sbox1(eng, input[88:96], vect_a, bit, ancilla[152:190])
    input[80:88] = Sbox2(eng, input[80:88], vect_b, bit, ancilla[190:228])
    input[72:80] = Sbox1Inv(eng, input[72:80], vect_a, bit, ancilla[228:266])
    input[64:72] = Sbox2Inv(eng, input[64:72], vect_b, bit, ancilla[266:304])

    input[56:64] = Sbox1(eng, input[56:64], vect_a, bit, ancilla[304:342])
    input[48:56] = Sbox2(eng, input[48:56], vect_b, bit, ancilla[342:380])
    input[40:48] = Sbox1Inv(eng, input[40:48], vect_a, bit, ancilla[380:418])
    input[32:40] = Sbox2Inv(eng, input[32:40], vect_b, bit, ancilla[418:456])

    input[24:32] = Sbox1(eng, input[24:32], vect_a, bit, ancilla[456:494])
    input[16:24] = Sbox2(eng, input[16:24], vect_b, bit, ancilla[494:532])
    input[8:16] = Sbox1Inv(eng, input[8:16], vect_a, bit, ancilla[532:570])
    input[0:8] = Sbox2Inv(eng, input[0:8], vect_b, bit, ancilla[570:])

    ## Diffusion layer
    input = Diffusion(eng, input)

    return input


def RoundEven(eng, input, key, vect_a, vect_b, bit, ancilla):
    ## input ^ key
    for i in range(128):
        CNOT | (key[i], input[i])

    ## S-box type2
    input[120:] = Sbox1Inv(eng, input[120:], vect_a, bit, ancilla[:38])
    input[112:120] = Sbox2Inv(eng, input[112:120], vect_b, bit, ancilla[38:76])
    input[104:112] = Sbox1(eng, input[104:112], vect_a, bit, ancilla[76:114])
    input[96:104] = Sbox2(eng, input[96:104], vect_b, bit, ancilla[114:152])

    input[88:96] = Sbox1Inv(eng, input[88:96], vect_a, bit, ancilla[152:190])
    input[80:88] = Sbox2Inv(eng, input[80:88], vect_b, bit, ancilla[190:228])
    input[72:80] = Sbox1(eng, input[72:80], vect_a, bit, ancilla[228:266])
    input[64:72] = Sbox2(eng, input[64:72], vect_b, bit, ancilla[266:304])

    input[56:64] = Sbox1Inv(eng, input[56:64], vect_a, bit, ancilla[304:342])
    input[48:56] = Sbox2Inv(eng, input[48:56], vect_b, bit, ancilla[342:380])
    input[40:48] = Sbox1(eng, input[40:48], vect_a, bit, ancilla[380:418])
    input[32:40] = Sbox2(eng, input[32:40], vect_b, bit, ancilla[418:456])

    input[24:32] = Sbox1Inv(eng, input[24:32], vect_a, bit, ancilla[456:494])
    input[16:24] = Sbox2Inv(eng, input[16:24], vect_b, bit, ancilla[494:532])
    input[8:16] = Sbox1(eng, input[8:16], vect_a, bit, ancilla[532:570])
    input[0:8] = Sbox2(eng, input[0:8], vect_b, bit, ancilla[570:])

    ## Diffusion layer
    input = Diffusion(eng, input)

    return input


def RoundFinal(eng, input, key, vect_a, vect_b, bit, ancilla):
    ## input ^ key
    for i in range(128):
        CNOT | (key[i], input[i])

    ## S-box type2
    input[120:] = Sbox1Inv(eng, input[120:], vect_a, bit, ancilla[:38])
    input[112:120] = Sbox2Inv(eng, input[112:120], vect_b, bit, ancilla[38:76])
    input[104:112] = Sbox1(eng, input[104:112], vect_a, bit, ancilla[76:114])
    input[96:104] = Sbox2(eng, input[96:104], vect_b, bit, ancilla[114:152])

    input[88:96] = Sbox1Inv(eng, input[88:96], vect_a, bit, ancilla[152:190])
    input[80:88] = Sbox2Inv(eng, input[80:88], vect_b, bit, ancilla[190:228])
    input[72:80] = Sbox1(eng, input[72:80], vect_a, bit, ancilla[228:266])
    input[64:72] = Sbox2(eng, input[64:72], vect_b, bit, ancilla[266:304])

    input[56:64] = Sbox1Inv(eng, input[56:64], vect_a, bit, ancilla[304:342])
    input[48:56] = Sbox2Inv(eng, input[48:56], vect_b, bit, ancilla[342:380])
    input[40:48] = Sbox1(eng, input[40:48], vect_a, bit, ancilla[380:418])
    input[32:40] = Sbox2(eng, input[32:40], vect_b, bit, ancilla[418:456])

    input[24:32] = Sbox1Inv(eng, input[24:32], vect_a, bit, ancilla[456:494])
    input[16:24] = Sbox2Inv(eng, input[16:24], vect_b, bit, ancilla[494:532])
    input[8:16] = Sbox1(eng, input[8:16], vect_a, bit, ancilla[532:570])
    input[0:8] = Sbox2(eng, input[0:8], vect_b, bit, ancilla[570:])

    return input


def Diffusion(eng, a):
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


def CNOT8bit(a, an, b, bn):
    for i in range(8):
        CNOT | (a[127 - an * 8 - i], b[127 - bn * 8 - i])


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
    # 기존의 A
    P = [1, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 1,
         0, 0, 0, 0, 0, 1, 0, 0]

    L = [1, 0, 0, 0, 0, 0, 0, 0,
         1, 1, 0, 0, 0, 0, 0, 0,
         1, 1, 1, 0, 0, 0, 0, 0,
         1, 1, 1, 1, 0, 0, 0, 0,
         1, 1, 1, 1, 1, 0, 0, 0,
         0, 0, 0, 1, 1, 1, 0, 0,
         0, 1, 1, 1, 0, 0, 1, 0,
         0, 0, 1, 1, 1, 0, 0, 1]

    U = [1, 0, 0, 0, 1, 1, 1, 1,
         0, 1, 0, 0, 1, 0, 0, 0,
         0, 0, 1, 0, 0, 1, 0, 0,
         0, 0, 0, 1, 0, 0, 1, 0,
         0, 0, 0, 0, 1, 0, 0, 1,
         0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 1]

    Apply_U(eng, inv, n, U)
    Apply_L(eng, inv, n, L)

    out = []
    for i in range(5):
        out.append(inv[i])
    out.append(inv[6])
    out.append(inv[7])
    out.append(inv[5])

    Round_constant_XOR(eng, out, b, n) # out = out ^ b

    return out


def MatrixProductS1Inv(eng, a, b, n):
    Round_constant_XOR(eng, a, b, n) # a = a ^ b

    # 기존의 A^-1
    P = [0, 0, 1, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 1]

    L = [1, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         1, 0, 1, 1, 0, 0, 0, 0,
         0, 1, 0, 1, 1, 0, 0, 0,
         0, 0, 1, 0, 1, 1, 0, 0,
         1, 0, 0, 0, 0, 1, 1, 0,
         0, 1, 0, 0, 0, 0, 1, 1]

    U = [1, 0, 0, 1, 0, 0, 1, 0,
         0, 1, 0, 0, 1, 0, 0, 1,
         0, 0, 1, 0, 0, 1, 0, 1,
         0, 0, 0, 1, 0, 0, 1, 1,
         0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 1]

    Apply_U(eng, a, n, U)
    Apply_L(eng, a, n, L)

    out = []
    out.append(a[2])
    out.append(a[0])
    out.append(a[1])

    for i in range(3, 8):
        out.append(a[i])

    return out


def MatrixProductS2(eng, inv, b, n):
    # 기존의 B
    P = [0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 1]

    L = [1, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         1, 0, 1, 1, 0, 0, 0, 0,
         1, 0, 1, 1, 1, 0, 0, 0,
         0, 1, 0, 1, 0, 1, 0, 0,
         0, 1, 1, 0, 1, 1, 1, 0,
         1, 0, 0, 1, 1, 1, 1, 1]

    U = [1, 1, 1, 0, 1, 1, 0, 1,
         0, 1, 0, 1, 0, 1, 1, 1,
         0, 0, 1, 1, 1, 1, 1, 1,
         0, 0, 0, 1, 0, 0, 0, 1,
         0, 0, 0, 0, 1, 1, 0, 1,
         0, 0, 0, 0, 0, 1, 0, 1,
         0, 0, 0, 0, 0, 0, 1, 1,
         0, 0, 0, 0, 0, 0, 0, 1]

    Apply_U(eng, inv, n, U)
    Apply_L(eng, inv, n, L)

    out = []
    out.append(inv[1])
    out.append(inv[2])
    out.append(inv[0])
    out.append(inv[3])

    out.append(inv[5])
    out.append(inv[4])
    out.append(inv[6])
    out.append(inv[7])

    Round_constant_XOR(eng, out, b, n) # out = out ^ b
    return out


def MatrixProductS2Inv(eng, a, b, n):
    Round_constant_XOR(eng, a, b, n) # a = a ^ b

    P = [0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0,
         1, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 1,
         0, 0, 0, 0, 0, 0, 1, 0]

    L = [1, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0,
         1, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 0, 0,
         1, 1, 0, 1, 0, 0, 1, 0,
         1, 1, 1, 1, 1, 1, 0, 1]

    U = [1, 1, 1, 0, 0, 0, 1, 1,
         0, 1, 1, 0, 1, 0, 1, 1,
         0, 0, 1, 0, 0, 1, 1, 0,
         0, 0, 0, 1, 1, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 1,
         0, 0, 0, 0, 0, 1, 0, 1,
         0, 0, 0, 0, 0, 0, 1, 1,
         0, 0, 0, 0, 0, 0, 0, 1]

    Apply_U(eng, a, n, U)
    Apply_L(eng, a, n, L)

    out = []
    out.append(a[3])
    out.append(a[2])
    out.append(a[5])
    out.append(a[0])

    out.append(a[4])
    out.append(a[1])
    out.append(a[7])
    out.append(a[6])

    return out


def Apply_U(eng, key, n, U):
    for i in range(n - 1):
        for j in range(n - 1 - i):
            if (U[(i * n) + 1 + i + j] == 1):
                CNOT | (key[1 + i + j], key[i])


def Apply_L(eng, key, n, L):
    for i in range(n - 1):
        for j in range(n - 1, 0 + i, -1):
            if (L[n * (n - 1 - i) + (n - 1 - j)] == 1):
                CNOT | (key[n - 1 - j], key[n - 1 - i])

def recursive_karatsuba(eng, a, b, n, count, ancilla):  # n=8
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
    c_b, count, ancilla = recursive_karatsuba(eng, a[r_low:n], b[r_low:n], n // 2, count, ancilla)  # 2 qubits     # 3~5
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
    P = [1, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 1]

    L = [1, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 0,  # (x4, x6) (x5, x6)
         0, 0, 0, 0, 0, 0, 1, 1]  # (x6, x7)

    U = [1, 0, 0, 0, 1, 0, 1, 0,  # (x4, x0) (x6, x0)
         0, 1, 0, 0, 0, 1, 0, 0,  # (x5, x1)
         0, 0, 1, 0, 1, 0, 0, 1,  # (x4, x2) (x7, x2)
         0, 0, 0, 1, 0, 1, 0, 0,  # (x5, x3)
         0, 0, 0, 0, 1, 0, 1, 1,  # (x6, x4) (x7, x4)
         0, 0, 0, 0, 0, 1, 1, 0,  # (x6, x5)
         0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 1]

    Apply_U(eng, vector, n, U)
    Apply_L(eng, vector, n, L)

    out = []
    out.append(vector[0])
    out.append(vector[4])
    out.append(vector[1])
    out.append(vector[6])

    out.append(vector[2])
    out.append(vector[5])
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

MK = 0x00112233445566778899aabbccddeeff
P = 0x11111111aaaaaaaa11111111bbbbbbbb
NCT = 1
resource_check = 0
AND_check = 0

# 결과값 확인
# classic = ClassicalSimulator()
# eng = MainEngine(classic)
# print("#### Ciphertext ####")
# ARIA128(eng)
# eng.flush()
# print('\n')

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

