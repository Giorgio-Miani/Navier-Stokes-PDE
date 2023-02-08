#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import sys

convergence_data = pd.read_csv(sys.argv[1], sep = ",")

plt.rcParams.update({"font.size": 11})

lista = []
for i in range(65):
    lista.append(i)

plt.plot(lista,
         convergence_data.Cl,
         label = 'Cu')

plt.xlabel("T")
plt.ylabel("Cd")
plt.legend()

plt.savefig("convergence.pdf")