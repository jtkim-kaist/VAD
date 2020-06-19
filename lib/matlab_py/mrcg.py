import numpy as np
import numpy.matlib
try:
    from scipy.fftpack import fft, ifft
except ImportError:
    from numpy.fft import fft, ifft
from scipy.signal import lfilter
import scipy.io as sio
import time
from scipy import signal
import os

epsc = 0.000001


def mrcg_features( sig, sampFreq = 16000):
    # This function computes MRCG features
    beta = 1000 / np.sqrt(sum(map(lambda x:x*x,sig)) / len(sig))
    sig = sig*beta
    sig = sig.reshape(len(sig), 1)
    t0= time.clock()
    g = gammatone(sig, 64, sampFreq) #Gammatone filterbank responses
    t1 = time.clock()
    cochlea1 = np.log10(cochleagram(g, int(sampFreq * 0.025), int(sampFreq * 0.010)))
    t2=time.clock()
    cochlea2 = np.log10(cochleagram(g, int(sampFreq * 0.200), int(sampFreq * 0.010)))
    t3 = time.clock()
    print('gamma total')
    print(t1-t0)
    print('coch1')
    print(t2-t1)
    print('coch2')
    print(t3-t2)
    cochlea1 = cochlea1[:,:]
    cochlea2 = cochlea2[:,:]
    t4 = time.clock()
    cochlea3 = get_avg(cochlea1, 5, 5)
    cochlea4 = get_avg(cochlea1, 11, 11)
    t5 = time.clock()
    print('get avg')
    print(t5-t4)
    all_cochleas = np.concatenate([cochlea1,cochlea2,cochlea3,cochlea4],0)

    del0 = deltas(all_cochleas)
    ddel = deltas(deltas(all_cochleas, 5), 5)

    ouotput = np.concatenate((all_cochleas, del0, ddel), 0)

    return ouotput

def gammatone(insig, numChan=128, fs = 16000):
    # Produce an array of filtered responses from a Gammatone filterbank.
    # The first variable is required.
    # numChan: number of filter channels.
    # fRange: frequency range.
    # fs: sampling frequency.
    # Written by ZZ Jin, adapted by DLW in Jan'07 and JF Woodruff in Nov'08
    # default number of filter channels in filterbank is 128
    # default sampling frequency is 16000
    fRange = [50, 8000] # default frequency range in Hz
    filterOrder = 4 # filter order
    gL = 2048 # gammatone filter length or 128 ms for 16 kHz sampling rate
    sigLength = len(insig) # input signal length
    phase = np.zeros([numChan, 1]) # initial phases
    erb_b = hz2erb(fRange) # upper and lower bound of ERB

    erb_b_diff = (erb_b[1]-erb_b[0])/(numChan-1)
    erb = np.arange(erb_b[0], erb_b[1]+epsc, erb_b_diff) # ERB segment
    cf = erb2hz(erb) # center frequency array indexed by channel
    b = [1.019 * 24.7 * (4.37 * x / 1000 + 1) for x in cf] # rate of decay or bandwidth

    #  Generating gammatone impulse responses with middle-ear gain normalization
    gt = np.zeros([numChan, gL]) # Initialization
    tmp_t = np.arange(1,gL+1)/fs
    for i in range(numChan):
        gain = 10**((loudness(cf[i])-60)/20)/3*(2 * np.pi * b[i] / fs)**4 # loudness-based gain adjustments
        tmp_temp = [gain*(fs**3)*x**(filterOrder - 1)*np.exp(-2 * np.pi * b[i] * x)*np.cos(2 * np.pi * cf[i] * x + phase[i]) for x in tmp_t]
        tmp_temp2 = np.reshape(tmp_temp, [1, gL])

        gt[i, :] = tmp_temp2

    sig = np.reshape(insig,[sigLength,1]) # convert input to column vector
    gt2 = np.transpose(gt)
    resig = np.matlib.repmat(sig,1,numChan)
    t0 = time.clock()
    r = np.transpose(fftfilt(gt2,resig,numChan)) # gammatone filtering using FFTFILT
    t1 = time.clock()
    print('fftfilter')
    print(t1-t0)
    return r

def hz2erb(hz):
    # Convert ERB-rate scale to normal frequency scale.
    # Units are number of ERBs and number of Hz.
    # ERB stands for Equivalent Rectangular Bandwidth.
    # Written by ZZ Jin, and adapted by DLW in Jan'07
    erb1 = 0.00437
    erb2 = np.multiply(erb1,hz)
    erb3 = np.subtract(erb2,-1)
    erb4 = np.log10(erb3)
    erb = 21.4 *erb4
    return erb

def erb2hz(erb):
    # Convert normal frequency scale in hz to ERB-rate scale.
    # Units are number of Hz and number of ERBs.
    # ERB stands for Equivalent Rectangular Bandwidth.
    # Written by ZZ Jin, and adapted by DLW in Jan'07
    hz = [(10**(x/21.4)-1)/(0.00437) for x in erb]
    return hz

def loudness(freq):
    # Compute loudness level in Phons on the basis of equal-loudness functions.
    # It accounts a middle ear effect and is used for frequency-dependent gain adjustments.
    # This function uses linear interpolation of a lookup table to compute the loudness level,
    # in phons, of a pure tone of frequency freq using the reference curve for sound
    # pressure level dB. The equation is taken from section 4 of BS3383.
    # Written by ZZ Jin, and adapted by DLW in Jan'07
    dB=60
    if os.path.exists('./lib/matlab_py/f_af_bf_cf.mat'):
        fmat = sio.loadmat('./lib/matlab_py/f_af_bf_cf.mat')
    else:
        fmat = sio.loadmat('.library/VAD/lib/matlab_py/f_af_bf_cf.mat')
    # Stores parameters of equal-loudness functions from BS3383,"Normal equal-loudness level
    # contours for pure tones under free-field listening conditions", table 1.
    # f (or ff) is the tone frequency, af and bf are frequency-dependent coefficients, and
    # tf is the threshold sound pressure level of the tone, in dBs
    af = fmat['af'][0]
    bf = fmat['bf'][0]
    cf = fmat['cf'][0]
    ff = fmat['ff'][0]
    i = 0
    while ff[i] < freq:
        i = i + 1

    afy = af[i - 1] + (freq - ff[i - 1]) * (af[i] - af[i - 1]) / (ff[i] - ff[i - 1])
    bfy = bf[i - 1] + (freq - ff[i - 1]) * (bf[i] - bf[i - 1]) / (ff[i] - ff[i - 1])
    cfy = cf[i - 1] + (freq - ff[i - 1]) * (cf[i] - cf[i - 1]) / (ff[i] - ff[i - 1])
    loud = 4.2 + afy * (dB - cfy) / (1 + bfy * (dB - cfy))
    return loud


def nextpow2(x):
    """Return the first integer N such that 2**N >= abs(x)"""
    return np.ceil(np.log2(abs(x)))


def fftfilt(b,x,nfft):
    fftflops = [18, 59, 138, 303, 660, 1441, 3150, 6875, 14952, 32373, 69762,
                149647, 319644, 680105, 1441974, 3047619, 6422736, 13500637, 28311786,
                59244791, 59244791*2.09]
    nb, _ = np.shape(b)
    nx, mx = np.shape(x)
    n_min = 0
    while 2**n_min < nb-1:
        n_min = n_min+1
    n_temp = np.arange(n_min, 21 + epsc, 1)
    n = np.power(2,n_temp)
    fftflops = fftflops[n_min-1:21]
    L = np.subtract(n,nb-1)
    temp_ind0 = np.ceil(np.divide(nx,L))
    temp_ind = np.multiply(temp_ind0,fftflops)
    temp_ind = np.array(temp_ind)
    ind = np.argmin(temp_ind)
    nfft=int(n[ind])
    L=int(L[ind])
    b_tr = np.transpose(b)
    B_tr = fft(b_tr,nfft)
    B = np.transpose(B_tr)
    y = np.zeros([nx, mx])
    istart = 0
    while istart < nx :
        iend = min(istart+L,nx)
        if (iend - istart) == 1 :
            X = x[0][0]*np.ones([nx,mx])
        else :
            xtr = np.transpose(x[istart:iend][:])
            Xtr = fft(xtr,nfft)
            X = np.transpose(Xtr)
        temp_Y = np.transpose(np.multiply(B,X))
        Ytr = ifft(temp_Y,nfft)
        Y = np.transpose(Ytr)
        yend = np.min([nx, istart + nfft])
        y[istart:yend][:] = y[istart:yend][:] + np.real(Y[0:yend-istart][:])

        istart = istart + L
    return y


def cochleagram(r, winLength = 320, winShift=160):
    # Generate a cochleagram from responses of a Gammatone filterbank.
    # It gives the log energy of T-F units
    # The first variable is required.
    # winLength: window (frame) length in samples
    # winShift : frame shift (default is half frame)
    # Written by ZZ Jin, and adapted by DLW in Jan'07
    # default window length in sample points which is 20 ms for 16 KHz sampling frequency
    numChan, sigLength = np.shape(r) # number of channels and input signal length
    M = np.floor(sigLength / winShift) # number of time frames
    a = np.zeros([numChan, int(M)])
    rs = np.square(r)
    rsl = np.concatenate((np.zeros([numChan,winLength-winShift]),rs),1)
    # calculate energy for each frame in each channel
    for m in range(int(M)):
        temp = rsl[:,m*winShift : m*winShift+winLength]
        a[:, m] = np.sum(temp,1)

    return a


def get_avg( m , v_span, h_span):
    # This function produces a smoothed version of cochleagram
    fil_size = (2 * v_span + 1) * (2 * h_span + 1)
    meanfil = np.ones([1+2*h_span,1+2*h_span])
    meanfil = np.divide(meanfil,fil_size)

    out = signal.convolve2d(m, meanfil, boundary='fill', fillvalue=0, mode='same')
    return out


def deltas(x, w=9) :
    # D = deltas(X,W)  Calculate the deltas (derivatives) of a sequence
    # Use a W-point window (W odd, default 9) to calculate deltas using a
    # simple linear slope.  This mirrors the delta calculation performed
    # in feacalc etc.  Each row of X is filtered separately.
    # 2003-06-30 dpwe@ee.columbia.edu

    nr,nc = np.shape(x)
    if nc ==0 :
        d= x # empty vector passed in; return empty vector
    else :
        hlen = int(np.floor(w / 2))
        win=np.arange(hlen, int(-(hlen+1)), -1)
        temp = x[:, 0]
        fx = np.matlib.repmat(temp.reshape([-1,1]), 1, int(hlen))
        temp = x[:, nc-1]
        ex = np.matlib.repmat(temp.reshape([-1,1]), 1, int(hlen))
        xx = np.concatenate((fx, x, ex),1)
        d = lfilter(win, 1, xx, 1)
        d = d[:,2*hlen:nc+2*hlen]

    return d
