function spc=hz2spc(hz,fs,n)
spc=round(2*hz/fs*(n-1));