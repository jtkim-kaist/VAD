function featmvn = mvn(feat)
% Mean and Variance Normalization

feat_mean = mean(feat,2);
feat_std = sqrt(var(feat,1,2));
featmvn = (feat - repmat(feat_mean,1,size(feat,2))) ./ repmat(feat_std,1,size(feat,2));
