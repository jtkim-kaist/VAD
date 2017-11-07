function erb=hz2erb(hz)
% Convert normal frequency scale in hz to ERB-rate scale.
% Units are number of Hz and number of ERBs.
% ERB stands for Equivalent Rectangular Bandwidth.
% Written by ZZ Jin, and adapted by DLW in Jan'07

erb=21.4*log10(4.37e-3*hz+1);

