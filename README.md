# VAD_Toolkit
This toolkit provides the voice activity detection (VAD) code and our recorded dataset.
## Introduction

![alt tag](https://user-images.githubusercontent.com/24668469/32532813-2b9c59aa-c490-11e7-8a30-a39de5aedc98.jpg)

VAD toolkit in this project was used in the paper: 

Juntae Kim, Minsoo Hahan: "Voice Activity Detection Based on the Adaptive Context Attention Model", submitted paper, 2017. 
This paper will be provided as soon as it is accepted.

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

![alt tag](https://user-images.githubusercontent.com/24668469/32533149-5526a77e-c492-11e7-909f-a7c7983d9dd4.jpg)

## Recorded Dataset

[Link](http://sail.ipdisk.co.kr:80/publist/VOL1/Database/VAD_DB/Recorded_data.zip)

## Trouble Shooting
## References
