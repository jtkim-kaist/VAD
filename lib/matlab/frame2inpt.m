function [ raw_label ] = frame2inpt(label, win_len, win_step)

num_frame = length(label);

total_len = (num_frame-1)*win_step + win_len;

raw_label = zeros(1, total_len);
% temp_label = zeros(1,win_len);
start_indx = 0;
i = 1;
while(1)
    if(start_indx+win_len>total_len)
        break;
    end
    if i==1
        raw_label(start_indx+1 : start_indx + win_len) = label(i);
    else
        temp_label(:) = label(i);
        raw_label(start_indx+1 : start_indx + win_len) = (raw_label(start_indx+1 : start_indx + win_len) + temp_label)/2;
%         temp_label = zeros(1,total_len);
    end
    i = i + 1;
    start_indx = start_indx + win_step;
end
    raw_label = movmean(raw_label, 20);
% raw_label = raw_label>= 1;

end