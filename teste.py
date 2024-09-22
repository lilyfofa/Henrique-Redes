import numpy as np

x = [2, 1.5, 0.9, 0, 0, 0, 0, 0, 0]

rmse = np.sqrt(np.mean(np.array([(x[i]-1)**2 for i in range(3)])))

print(rmse)