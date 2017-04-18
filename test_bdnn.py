import numpy as np
import utils

a = np.arange(12) + 1
a = a.reshape((-1, 3))

bb = utils.bdnn_transform(a, 2, 1)

for i in range(bb.shape[0]):
    print(bb[i, :])




