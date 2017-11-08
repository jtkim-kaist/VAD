# Voice Activity Detection Toolkit
This toolkit provides the voice activity detection (VAD) code and our recorded dataset.
## Introduction

![alt tag](https://user-images.githubusercontent.com/24668469/32532813-2b9c59aa-c490-11e7-8a30-a39de5aedc98.jpg)

VAD toolkit in this project was used in the paper: 

Juntae Kim, Minsoo Hahan: "Voice Activity Detection Based on the Adaptive Context Attention Model", submitted paper, 2017.

This paper will be provided as soon as it is accepted. If you want to use this toolkit before publishing the paper, please cite this
repositoy like: 

Juntae Kim: "VAD_Toolkit", GitHub repository, [Online] Available: https://github.com/jtkim-kaist/VAD_Toolkit, 2017.

VAD in this toolkit follows the procedure as below:

#### Acoustic feature extraction

In this toolkit, we use the The multi-resolution cochleagram (MRCG) [1] for the acoustic feature implemented by matlab.
Note that MRCG extraction time is relatively long compared to the classifier.
#### Classifier

This toolkit supports 4 types of MRCG based classifer implemented by python with tensorflow as follows:
1. Adaptive context attention model (ACAM)
2. Boosted deep neural network (bDNN) [1]
3. Deep neural network (DNN) [1] 
4. Long short term memroy recurrent neural network (LSTM)

## Prerequisites

- Python 3

- Tensorflow 1.1

- Matlab 2017b
## Example

The example matlab script is `main.m`. Just run it on the matlab.
The result will be like following figure. 

Note: To apply this toolkit to other speech data, the speech data should be sampled with 16kHz sampling frequency.

![alt tag](https://user-images.githubusercontent.com/24668469/32533149-5526a77e-c492-11e7-909f-a7c7983d9dd4.jpg)

## Recorded Dataset
Our recored dataset is available: 
[Download](http://sail.ipdisk.co.kr:80/publist/VOL1/Database/VAD_DB/Recorded_data.zip)


#### Specification
- Environments

>Bus stop, construction site, park, and room.

- Recording device

>A smart phone (Samsung Galaxy 8)

At each environment, conversational speech by two Korean male speakers was recorded. The ground truth labels are manually annotated. Because the recording was carried out in the real world, unexpected noises are included to the dataset such as the crying of baby, the chirping of insects, mouse click sound, and etc. The details of dataset is described in the following table:


|               | Bus stop      | Cons. site    | Park          | Room          | Overall       |
| :------------ | :-----------: | :-----------: | :-----------: | :-----------: | :-----------: |
| Dur. (min)    | 30.02         | 30.03         | 30.07         | 30.05         | 120.17        |
| Avg. SNR (dB) | 5.61          | 2.05          | 5.71          | 18.26         | 7.91          |
| % of speech   | 40.12         | 26.71         | 26.85         | 30.44         | 31.03         |
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
[1] Zhang, Xiao-Lei, and DeLiang Wang. “Boosting contextual information for deep neural network based voice activity detection,” IEEE Trans. Audio, Speech, Lang. Process., vol. 24, no. 2, pp. 252-264, 2016.
