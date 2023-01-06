function [sf,r]  =csp(X1,X2,n)
% X1:EEG in condition 1 with format chan*time*trial1
% X2:EEG in condition 1 with format chan*time*trial2
% n: number of spatial filters

for iTrail = 1:size(X1, 1)
    for iCh = 1:size(X1, 2)
        X1Trans(iCh, :, iTrail) = X1(iTrail, iCh, :);
    end
end

for iTrail = 1:size(X2, 1)
    for iCh = 1:size(X2, 2)
        X2Trans(iCh, :, iTrail) = X2(iTrail, iCh, :);
    end
end

X1 = X1Trans;
X2 = X2Trans;

if nargin<3
    n = 1;
end
[chan1,len1,trial1]  = size(X1);
[chan2,len2,trial2]  =size(X2);
if (chan1 ~= chan2)
    error('channel of X1 do not equal to that of X2!');
end

X1_mean = mean(X1,3);
X2_mean = mean(X2,3);

covX1 = X1(:,:)-repmat(X1_mean,1,trial1);
covX1 = covX1*covX1'/trial1;

covX2 = X2(:,:)-repmat(X2_mean,1,trial2);
covX2 = covX2*covX2'/trial2;

[V,D] = eig(covX1 ,covX1+ covX2);

r = abs(diag(D));
[Y,I] = sort(r,'descend');
V = V(:,I);
r = Y;
r = [r(1:n),r(chan1-n+1:chan1)];
sf =[V(:,1:n),V(:,chan1-n+1:chan1)];
