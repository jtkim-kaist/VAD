function [ Detect ] = Truelabel2Trueframe( TrueLabel_bin,wsize,wstep )


iidx=1;
Frame_iidx=1;
Frame_len=Frame_Length( TrueLabel_bin,wstep,wsize );
Detect=zeros(Frame_len,1);

while(1)
    
    %% Windowing
    if(iidx+wsize<length(TrueLabel_bin))
        TrueLabel_frame=TrueLabel_bin(iidx:iidx+wsize-1)*10; % Because the size of TrueLabel_bin is 0.1 in TrueLabel function
    else
        TrueLabel_frame=TrueLabel_bin(iidx:end)*10; % Because the size of TrueLabel_bin is 0.1 in TrueLabel function
    end
    
    if(sum(TrueLabel_frame)>=wsize/2)
        TrueLabel_frame=1;
    else
        TrueLabel_frame=0;
    end
    if(Frame_iidx>length(Detect))
        break;
    end
    Detect(Frame_iidx)=TrueLabel_frame;
    iidx=iidx+wstep;
    Frame_iidx=Frame_iidx+1;
    if(iidx>length(TrueLabel_bin))
        break;
    end
    
end

end

