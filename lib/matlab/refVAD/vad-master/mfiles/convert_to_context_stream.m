function feat = convert_to_context_stream( x , K, order)

if ~exist('K', 'var')
   K=30; % 300 ms
end
if ~exist('order', 'var')
   order=4;
end
[M,Mstatic]=vec2featmat_d1(K,order);

% select static + delta's
nbp=length(x)-K+1;
ind1=(1:nbp)-1;
ind2=(1:K)';
x_aug=x(ind1(ones(K,1),:)+ind2(:,ones(1,nbp)));

feat=Mstatic*x_aug;
feat=[feat(:,nbp*ones(1,K/2)) feat feat(:,nbp*ones(1,K/2))];
feat=feat(:,1:length(x));
