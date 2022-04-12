import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import pandas as pd

u = 0
w = 1

if u == 0:
    if w > 0:
        phi = np.pi / 2
    else:
        phi = -np.pi / 2
else:
    phi = np.arctan2(w,u)

if u < 0:
    phi = phi + np.pi

print(phi)