#!/bin/sh

curdir=`pwd`

rm -rf $curdir/data/feat/*

rm -rf $curdir/logs/DNN/*
rm -rf $curdir/logs/bDNN/*
rm -rf $curdir/logs/LSTM/*
rm -rf $curdir/logs/ACAM/*

rm -rf $curdir/norm_data/*.mat
rm -rf $curdir/result/*

rm -rf $curdir/saved_model/graph/DNN/*
rm -rf $curdir/saved_model/graph/bDNN/*
rm -rf $curdir/saved_model/graph/LSTM/*
rm -rf $curdir/saved_model/graph/ACAM/*

rm -rf $curdir/sample_data/*



