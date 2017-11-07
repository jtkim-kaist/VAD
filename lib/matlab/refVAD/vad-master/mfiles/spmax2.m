function [y,arg]=spmax2(row,col,val1,val2,N)
% MATLAB version 
% [y,arg]=spmax2(row,col,val1,val2,N)
% is equivalent to [y,arg]=spmax(x)
% where x has N cols and is given by val1+val2 in position(row,col). If val1==-INF => entry disregarded
% -INF in y signify that column was empty.
y=-inf*ones(1,N);
arg=zeros(1,N);
val1=val1(:);
val2=val2(:);
interesting=find(val1~=-inf);
val=val1(interesting)+val2(interesting);
col=col(interesting);
row=row(interesting);
for i=1:N,
   sel=find(col==i);
   if ~isempty(sel),
      [y(i),arg(i)]=max(val(sel));
      arg(i)=row(sel(arg(i)));
   end
end
