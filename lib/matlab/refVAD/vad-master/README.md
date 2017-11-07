#Voice Activity Detection system

## Description
VAD system based on Deep Neural Networks (DNN) and feature fusion (Gammatone, Gabor, Long-term Spectral Variability and voicing).  
System was developed as part of the RATS (Robust Automated Transcription of Speech) program of DARPA.

##Instructions:

1. open Matlab
2. run the script **apply_vad(_path/to/audio_)**: 
3. a figure will appear that shows the original signal, and VAD labels, given a directory of audio wav files.
4. to control the accuracy (depending on how noisy the files are), you can play with the parameters p1 and p2  
5. additional info 

       Apply Voice Activity Detection to all files in a specified audio directory

       --IN--  
       audiodir: directory of audio files (WAV format)  
       p1: speech/non-speech threshold [default:0.1]  
       p2: speech region smoothing [default:20]  

       --OUT--  
       vadout: VAD labels at frame level (10 ms)  


##Reference citation:

Van Segbroeck, Maarten, Andreas Tsiartas, and Shrikanth Narayanan. _"A robust frontend for VAD: exploiting contextual, discriminative and spectral cues of human voice."_ INTERSPEECH. 2013.

**bibtex**

@inproceedings{van2013robust,
  title={A robust frontend for VAD: exploiting contextual, discriminative and spectral cues of human voice.},
  author={Van Segbroeck, Maarten and Tsiartas, Andreas and Narayanan, Shrikanth},
  booktitle={INTERSPEECH},
  pages={704--708},
  year={2013}
}
