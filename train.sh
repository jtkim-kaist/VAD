#!/bin/sh

curdir=`pwd`
pylibdir=`realpath $curdir/lib/python`
train=`realpath $pylibdir/train.py`
ckpt_update=`realpath $pylibdir/update_ckpt.py`
# train script options
# m 0 : DNN
# e : extract MRCG feature (1) or not (0)

python3 $train -m 0 -e 1 --train_step=100 --batch_size=256 --prj_dir=$curdir

a
# ckpt_update script options
# u : update checkpoint from trained model (1) or restore checkpoint to default (0)
# Note that when u==0, the normalization factor is also restored to default.
python3 $ckpt_update -u 1 --model=DNN --prj_dir=$curdir
