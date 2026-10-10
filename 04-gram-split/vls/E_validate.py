"""Validate the E data at Q = 6000 against D_numerics/data/A_Q6000_K12.npy (scales Y = 2^8..2^19)."""
import numpy as np
GS = "/Users/akashlevy/Documents/riemann/gram_split"
A = np.load(GS + "/vls/data/A_Q6000_E21.npy"); m = np.load(GS + "/vls/data/m_Q6000_E21.npy")
A0 = np.load(GS + "/D_numerics/data/A_Q6000_K12.npy"); m0 = np.load(GS + "/D_numerics/data/m_Q6000_K12.npy")
key = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(m[0], m[1]))}
ix = np.array([key[(int(a), int(b))] for a, b in zip(m0[0], m0[1])])
print("rows: E %d, K12 %d" % (A.shape[1], A0.shape[1]))
for i in range(A0.shape[0]):
    print("Y=2^%d: max|A_E - A_K12| = %.2e   (max|A| = %.2e)" % (8 + i, np.max(np.abs(A[i, ix] - A0[i])), np.max(np.abs(A0[i]))))
