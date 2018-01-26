#!/bin/sh

curdir=`pwd`
pylibdir=`realpath $curdir/lib/python`
train=`realpath $pylibdir/train.py`
# train script options
# m 0 : ACAM
# m 1 : bDNN
# m 2 : DNN
# m 3 : LSTM
# e : extract MRCG feature (1) or not (0)

python3 $train -m 0 -e 1 --prj_dir=$curdir

