import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd

filename = sys.argv[1]
df = pd.read_csv(filename)
df.insert(0, "p_logic", df["n_success"]/df["n_trials"])
df.insert(0, "N", 2*df["w"]*df["h"])

for N in df["N"].unique():
    df2 = df[df["N"] == N]
    label = "N=" + str(N)
    plt.plot(df2["p"],df2["p_logic"],label=label)

plt.xlabel("p")
plt.ylabel("p_logic")
plt.legend()
#plt.title("Depolarising Noise")
plt.show()

