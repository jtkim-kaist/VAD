import librosa
import numpy as np
import matplotlib.pyplot as plt

import lib.matlab_py.utils as utils

audio_dir = './data/example/clean_speech.wav'

mode = 0
th = 0.4

output_type = 1
is_default = 1

result = utils.vad_func(audio_dir, mode, th, output_type, is_default, off_on_length=20, on_off_length=20,
                        hang_before=20, hang_over=20)
s, audio_sr = librosa.load(audio_dir, sr=16000)

t_max = np.minimum(s.shape[0], result.shape[0])

t = np.arange(0, t_max)/audio_sr
s = np.divide(s, np.max(np.absolute(s)))

result = np.multiply(result, 0.3)
plt.plot(t, s[0:t_max], 'b')
plt.plot(t, result, 'g')
plt.show()
