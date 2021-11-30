import numpy as np

a = np.array([1,2,3])
b = np.zeros((3))
c = np.amax([a,b], axis=0)
print(c)
