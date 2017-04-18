import scipy.io as sio
import numpy as np

stft = sio.loadmat('nx_stft_003.mat')
stft = stft['s']
label = np.fromfile('./label_400.bin',dtype = np.float32)

print(np.shape(stft))
print(np.shape(label))
