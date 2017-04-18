# prediction


import numpy as np
import utils

batch_size = 20

idx = np.arange(batch_size) + 1
idx = np.reshape(idx, (batch_size, 1))


data = np.arange(batch_size*3) -1
data = np.reshape(data, (batch_size, 3))

label = np.arange(batch_size) + 10
label = label.reshape((batch_size,1))

w = 5
u = 2


trans_idx = utils.bdnn_transform(idx, w, u)
trans_data = utils.bdnn_transform(data, w, u)
trans_label = utils.bdnn_transform(label, w, u)

result = np.zeros((batch_size, 1))


final_data = trans_data[w:(batch_size-w), :]

numm  = np.arange(w, batch_size-w)

# for i in numm:
#     idx_temp = np.where((trans_idx-1) == i)
#     aa = trans_label[idx_temp]
#     aa = np.sum(aa)/aa.shape[0]
#     result[i] = aa

for i in range(batch_size):
    idx_temp = np.where((trans_idx-1) == i)
    aa = trans_label[idx_temp]
    aa = np.sum(aa)/aa.shape[0]
    result[i] = aa

ttt = np.trim_zeros(result)
print("done")
