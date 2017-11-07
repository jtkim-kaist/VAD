function [er, bad, out] = nntest(nn, x, y)
    if ~exist('y', 'var')
        y = zeros(size(x,1),nn.size(end));
    end
    nn.testing = 1;
    nn = nnff(nn, x, y);
    nn.testing = 0;
    
    [dum i] = max(nn.a{end},[],2);
    [dum, g] = max(y,[],2);
    bad = find(i ~= g);    
    er = numel(bad) / size(x, 1);
    out = nn.a{end};
end
