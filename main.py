import sys
sys.path.append("./lib/python")
#import librosa
import numpy as np
import matplotlib.pyplot as plt
import lib.mrcg.utils as util
import scipy.io.wavfile

audio_dir = './data/example/SNR103F3MIC021002_ch01.wav'

mode = 1
th= 0.4
output_type = 1
is_default = 1
result = util.vad_func(audio_dir, mode, th, output_type, is_default)
#s, audio_sr = librosa.load(audio_dir, sr=16000)
#result = util.vad_func(audio_dir, mode, th, output_type)
audio_sr, s = scipy.io.wavfile.read(audio_dir)
print(len(s)/audio_sr)
t_max = np.min((len(s),len(result[0])))
t = np.arange(0.,t_max/16000,1/16000)
s = np.divide(s,np.max(np.absolute(s)))
result = np.multiply(result,0.3)
plt.plot(t,s[0:t_max],'b')
plt.plot(t,result[0][0:t_max],'g')
plt.show()
