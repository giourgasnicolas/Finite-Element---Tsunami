import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def read_tsunami(filename):
    df = pd.read_csv(filename, sep=';')
    df.drop(df.columns[len(df.columns)-1], axis=1, inplace=True)
    df.drop("number", axis=1, inplace=True)
    df.drop("x", axis=1, inplace=True)
    df.drop("y", axis=1, inplace=True)

    E = df[df["name"] == 0]
    U = df[df["name"] == 1]
    V = df[df["name"] == 2]

    return (E["val"].to_numpy(), U["val"].to_numpy(), V["val"].to_numpy())

def extrapolate_rirchardson(FA, FB, i, j):
    E = ((2**j)*FB[0] - FA[0])/(2**j - 1)
    U = ((2**j)*FB[1] - FA[1])/(2**j - 1)
    V = ((2**j)*FB[2] - FA[2])/(2**j - 1)
    return (E, U, V)

def errors(Fij, Fref) :
    return (np.sum(np.abs(Fij[0] - Fref[0])),
            np.sum(np.abs(Fij[1] - Fref[1])),
            np.sum(np.abs(Fij[2] - Fref[2])))

def convergence(errA, errB) :
    conv_E = np.log2(errA[0] / errB[0])
    conv_U = np.log2(errA[1] / errB[1])
    conv_V = np.log2(errA[2] / errB[2])
    print("\n== Taux de convergence == ")
    print("==   pour E : {}".format(conv_E))
    print("==   pour U : {}".format(conv_U))
    print("==   pour V : {}".format(conv_V))
    return (conv_E, conv_U, conv_V)


F00 = read_tsunami("med_8.txt")
F10 = read_tsunami("med_4.txt")
F20 = read_tsunami("med_2.txt")
F30 = read_tsunami("med_1.txt")

F11 = extrapolate_rirchardson(F00, F10, 1, 1)
F21 = extrapolate_rirchardson(F10, F20, 2, 1)
F22 = extrapolate_rirchardson(F11, F21, 2, 2)
F31 = extrapolate_rirchardson(F20, F30, 3, 1)
F32 = extrapolate_rirchardson(F21, F31, 3, 2)
F33 = extrapolate_rirchardson(F22, F32, 3, 3)

ref = F33

err00 = errors(F00, ref)
err10 = errors(F10, ref)
err20 = errors(F20, ref)
err30 = errors(F30, ref)

C1 = convergence(err00, err10)
C2 = convergence(err10, err20)
C3 = convergence(err20, err30)

mean = 0
count = 0

for t in np.array([C1, C2, C3]).ravel() :
    mean += t
    count += 1

mean /= count

print("\nMoyenne des taux : {}\n".format(mean))
