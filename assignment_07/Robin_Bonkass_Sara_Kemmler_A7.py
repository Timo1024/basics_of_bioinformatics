from math import exp
import matplotlib.pyplot as plt
import numpy as np

def f(x, t):
    return x * exp(1) ** (-x * (1 - t))

x = np.linspace(0, 20, num=1000)

functions = {"red": 0.2, "blue": 0.3, "green": 0.4, "orange": 0.5, "black": 0.6}

plt.figure(figsize=(8, 6))

for i in functions:
    t = functions[i]
    text_label = "\u03B8 = " + str(t)
    plt.plot(x, f(x, t), color=i, linewidth=0.8, label=text_label)

plt.ylabel("expected number of contigs scaled in units G/L")
plt.xlabel("coverage c")
plt.legend()

plt.savefig("Robin_Bonkass_Sara_Kemmler_A7_plot.pdf", dpi=300, transparent=True)
