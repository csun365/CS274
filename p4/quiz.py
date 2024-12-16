import numpy as np
import matplotlib.pyplot as plt

with open("100.log", "r") as f:
    p_100 = [float(i) for i in f.readlines()]
with open("500.log", "r") as f:
    p_500 = [float(i) for i in f.readlines()]
with open("1000.log", "r") as f:
    p_1000 = [float(i) for i in f.readlines()]

print(np.mean(p_100), np.std(p_100))
print(np.mean(p_500), np.std(p_500))
print(np.mean(p_1000), np.std(p_1000))

plt.hist(p_100)
plt.show()

plt.hist(p_500)
plt.show()

plt.hist(p_1000)
plt.show()