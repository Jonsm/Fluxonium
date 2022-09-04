import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd

filename = sys.argv[1]
df = pd.read_csv(filename)
df.insert(0, "p_logic", df["n_success"]/df["n_trials"])

print(df)

coord_keys = ["l1_x","l1_y","l2_x","l2_y"]
for ls in np.unique(df[coord_keys].values,axis=0):
    df2 = df.copy()
    for i,l in enumerate(ls):
        df2 = df2[df2[coord_keys[i]]==l]
        
    label=f'l1=({ls[0]},{ls[1]}),l2=({ls[2]},{ls[3]})'
    plt.plot(df2["p"],df2["p_logic"],label=label)

plt.xlabel("p")
plt.ylabel("p_logic")
plt.legend()
plt.show()

