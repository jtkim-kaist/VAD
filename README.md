# Voice Activity Detection Toolkit
This toolkit provides the voice activity detection (VAD) code and our recorded dataset.
## Update
#### 2019-02-11

- Accepted for the presentation of this toolkit in ICASSP 2019!

#### 2018-12-11

- Post processing is updated. 

#### 2018-06-04

- Good news! we have uploaded speech enhancement toolkit based on deep neural network. This toolkit provides several useful things such as data generation script. You can find this toolkit in [here](https://github.com/jtkim-kaist/Speech-enhancement) 


#### 2018-04-09


- The test sciprt fully written by python has been uploaded in 'py' branch.


## Introduction

![alt tag](https://user-images.githubusercontent.com/24668469/32532813-2b9c59aa-c490-11e7-8a30-a39de5aedc98.jpg)

VAD toolkit in this project was used in the paper: 

J. Kim and M. Hahn, "Voice Activity Detection Using an Adaptive Context Attention Model," in *IEEE Signal Processing Letters*, vol. PP, no. 99, pp. 1-1.

URL: https://ieeexplore.ieee.org/document/8309294/

If our VAD toolkit supports your research, we are very appreciated if you cite this paper.

ACAM is based on the recurrent attention model (RAM) [1] and the implementation of RAM can be found in [jlindsey15](https://github.com/jlindsey15/RAM) and [jtkim-kaist](https://github.com/jtkim-kaist/ram_modified)'s repository.

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

- Matlab 2017b (will be depreciated)

## Example

The default model provided in this toolkit is the trained model using our dataset. The used dataset is described in our submitted paper.
The example matlab script is `main.m`. Just run it on the matlab. 
The result will be like following figure. 

Note: To apply this toolkit to other speech data, the speech data should be sampled with 16kHz sampling frequency.

![alt tag](https://user-images.githubusercontent.com/24668469/32533149-5526a77e-c492-11e7-909f-a7c7983d9dd4.jpg)

## Post processing

Many people want to the post-processing so I updated.

In py branch, you can see some parameters in utils.vad_func in main.py

Each parameter can handle following errors.



![alt tag](https://user-images.githubusercontent.com/24668469/49742392-cd778680-fcdb-11e8-96b9-a599a4f85f4f.PNG)


FEC: hang_before

MSC: off_on_length

OVER: hang_over

NDS: on_off_length

Note that there is NO optimal one. The optimal parameter set is according to the application.

Enjoy.

## Training
1. We attached the sample database to 'path/to/project/data/raw'. Please refer to the database for understanding the data format. 
2. The model specifications are described in 'path/to/project/configure'.
3. The training procedure has 2 steps: (i) MRCG extraction; (ii) Model training.

Note: Do not forget adding the path to this project in the matlab.

```
# train.sh
# train script options
# m 0 : ACAM
# m 1 : bDNN
# m 2 : DNN
# m 3 : LSTM
# e : extract MRCG feature (1) or not (0)

python3 $train -m 0 -e 1 --prj_dir=$curdir
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
1. Although MRCG show good performance but extraction time is somewhat long, therefore we will substitute it to other feature such as spectrogram.
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

## Acknowledgement

Jaeseok, Kim (KAIST) contributed to this project for changing matlab script to python.
