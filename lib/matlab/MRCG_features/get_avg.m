function out = get_avg( m , v_span, h_span)
%This function produces a smoothed version of cochleagram
[nr,nc] = size(m);
out = zeros(nr,nc);
fil_size = (2*v_span+1)*(2*h_span+1);

for i=1:nr
    row_begin = 1;
    row_end = nr;
    col_begin =1;
    col_end = nc;
    if (i - v_span)>=1
        row_begin = i - v_span;
    end    
    if (i + v_span)<=nr
        row_end = i + v_span;
    end
    
    for j=1:nc     
        if (j - h_span)>=1
            col_begin = j - h_span;
        end
        
        if (j + h_span)<=nc
            col_end = j + h_span;
        end
        
        tmp = m(row_begin:row_end,col_begin:col_end); 
        out(i,j) = sum(sum(tmp))/fil_size;
    end  
end

