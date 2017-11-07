function [ framelen ] = Frame_Length( x,overlap,nwind )
nx=length(x);
noverlap=nwind-overlap;
ncol = fix((nx-noverlap)/(nwind-noverlap));
colindex = 1 + (0:(ncol-1))*(nwind-noverlap);
framelen=length(colindex);
end

