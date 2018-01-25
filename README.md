# Voice Activity Detection Toolkit
This toolkit provides the voice activity detection (VAD) code and our recorded dataset.
## Introduction

![alt tag](https://user-images.githubusercontent.com/24668469/32532813-2b9c59aa-c490-11e7-8a30-a39de5aedc98.jpg)

VAD toolkit in this project was used in the paper: 

Juntae Kim, Minsoo Hahan: "Voice Activity Detection Based on the Adaptive Context Attention Model (ACAM)", submitted paper, 2017.

ACAM is based on the recurrent attention model (RAM) [1] and the implementation of RAM can be found in [jlindsey15](https://github.com/jlindsey15/RAM) and [jtkim-kaist](https://github.com/jtkim-kaist/ram_modified)'s repository.

The submitted paper will be provided as soon as it is accepted. If you want to use this toolkit before publishing the paper, please cite this repositoy like: 

Juntae Kim: "VAD_Toolkit", GitHub repository, [Online] Available: https://github.com/jtkim-kaist/VAD_Toolkit, 2017.

VAD in this toolkit follows the procedure as below:

#### Acoustic feature extraction

In this toolkit, we use the multi-resolution cochleagram (MRCG) [2] for the acoustic feature implemented by matlab.
Note that MRCG extraction time is relatively long compared to the classifier.
#### Classifier

This toolkit supports 4 types of MRCG based classifer implemented by python with tensorflow as follows:
1. Adaptive context attention model (ACAM)
2. Boosted deep neural network (bDNN) [2]
3. Deep neural network (DNN) [2] 
4. Long short term memory recurrent neural network (LSTM-RNN) [3]

## Prerequisites

- Python 3

- Tensorflow 1.1-3

- Matlab 2017b
## Example

The example matlab script is `main.m`. Just run it on the matlab.
The result will be like following figure. 

Note: To apply this toolkit to other speech data, the speech data should be sampled with 16kHz sampling frequency.

![alt tag](https://user-images.githubusercontent.com/24668469/32533149-5526a77e-c492-11e7-909f-a7c7983d9dd4.jpg)
## Training
We attached the sample database to 'path/to/project/data/raw'. Please refer to the database for understanding the data format. 
The training procedure has 2 steps: (i) MRCG extraction; (ii) Model training.

Note: Do not forget adding the path to this project in the matlab. Current version only supports DNN based training. We will update training script for other models.

```
# train.sh
# train script options
# m 0 : DNN
# e : extract MRCG feature (1) or not (0). 
# The MRCG extraction time is somewhat long so you can pass the feature extraction step if you already have MRCG feature.

python3 $train -m 0 -e 1 --train_step=100 --prj_dir=$curdir

# ckpt_update script options
# u : update checkpoint from trained model (1) or restore checkpoint to default (0).
# Note that when u==0, the normalization factor is also restored to default.
# After training you should update the model checkpoint with the normalization factor.

python3 $ckpt_update -u 1 --model=DNN --prj_dir=$curdir
```

## Recorded Dataset
Our recored dataset is freely available: 
[Download](http://sail.ipdisk.co.kr:80/publist/VOL1/Database/VAD_DB/Recorded_data.zip)


#### Specification
- Environments

>Bus stop, construction site, park, and room.

- Recording device

>A smart phone (Samsung Galaxy S8)

At each environment, conversational speech by two Korean male speakers was recorded. The ground truth labels are manually annotated. Because the recording was carried out in the real world, unexpected noises are included to the dataset such as the crying of baby, the chirping of insects, mouse click sound, and etc. The details of dataset is described in the following table:


|               | Bus stop      | Cons. site    | Park          | Room          | Overall       |
| :------------ | :-----------: | :-----------: | :-----------: | :-----------: | :-----------: |
| Dur. (min)    | 30.02         | 30.03         | 30.07         | 30.05         | 120.17        |
| Avg. SNR (dB) | 5.61          | 2.05          | 5.71          | 18.26         | 7.91          |
| % of speech   | 40.12         | 26.71         | 26.85         | 30.44         | 31.03         |
## TODO List
1. Freezing the graph for running the model fast.
2. Training script for bDNN, LSTM, ACAM --> will be updated until 2018-01-26
3. Although MRCG show good performance but extraction time is somewhat long, therefore we will substitute it to other feature such as spectrogram.
## Trouble Shooting
If you find any errors in the code, please contact to us.

E-mail: jtkim@kaist.ac.kr
## Copyright
Copyright (c) 2017 Speech and Audio Information Laboratory, KAIST, South Korea

License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
## References
[1] J. Ba, V. Minh, and K. Kavukcuoglu, “Multiple object recognition with visual attention,” arXiv preprint arXiv, 1412.7755, 2014.

[2] Zhang, Xiao-Lei, and DeLiang Wang. “Boosting contextual information for deep neural network based voice activity detection,” IEEE Trans. Audio, Speech, Lang. Process., vol. 24, no. 2, pp. 252-264, 2016.

[3] Zazo Candil, Ruben, et al. “Feature learning with raw-waveform CLDNNs for Voice Activity Detection.”, 2016.
