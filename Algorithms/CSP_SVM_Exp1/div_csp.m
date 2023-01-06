function [P,y]  =div_csp(X1,X2,n)
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

% X1_mean = mean(X1,3);
% X2_mean = mean(X2,3);

covX1 = zeros(size(X1,1),size(X1,1));
for iTrial = 1:trial1
    cov = X1(:,:,iTrial)*X1(:,:,iTrial)';
    allcovs{1}(:,:,iTrial) = cov/trace(cov);
end

covX2 = zeros(size(X2,1),size(X2,1));
for iTrial = 1:trial2
    cov = X2(:,:,iTrial)*X2(:,:,iTrial)';
    allcovs{2}(:,:,iTrial) = cov/trace(cov);
end

%%% configuration
opts.max_iter = 1000;					%% maximum number of iterations
opts.nreps = 1;							%% number of initializations
opts.lam = 0;								%% regularization parameter lambda
opts.beta = 0;								%% beta parameter (if 0 apply KL divergence, otherwise apply beta divergence)
opts.mode = 0;								%% if 0 substract regularization term (divCSP-WS, divCSP-BS, divCSP-AS), if 1 add regularization term (divCSP-M)
opts.deflation = 0;						%% if 0 then apply subspace algorithm, if 1 apply deflation method
opts.pca = 1;								%% find meaningful basis in extracted subspace by applying PCA
opts.csp_init = 1;						%% initialize first repetition with CSP solution
opts.sym = 1;								%% if 1 then apply symmetric divergence for regularization term, if 0 use the asymmetric divergence
opts.quiet = 1;							%% if 0 then show objective value for each step, if 1 then do not show
dim = 2*n; 									%% dimensionality of the CSP subspace (number of spatial filters)

%%% create epochs
cspTerm{1} = mean(allcovs{1},3);		%% average covariance matric for class 1
cspTerm{2} = mean(allcovs{2},3);		%% average covariance matric for class 2

% for i=1:size(allcovs{1},3)
% 	regTerm{1}(:,:,i) = allcovs{1}(:,:,i);
% 	regTerm{2}(:,:,i) = allcovs{2}(:,:,i);
% end
% opts.sym = 1;								%% use symmetric divergence
% opts.lam = 1;								%% only use the regularization term (sum of trial divergences)
% opts.mode = 1;								%% maximize regularization term instead of minimizing
% opts.beta = 0.1;							%% use beta divergence
% [P, results, y] = divcsp(cspTerm, regTerm, [], dim, opts);
% P = P';

%% for divCSP-WS the first term is the trial covariance matrix, the second term is the class average
for i=1:size(allcovs{1},3)
	regTerm{1}(:,:,i) = allcovs{1}(:,:,i);
	regTerm{2}(:,:,i) = cspTerm{1};
end
for i=1:size(allcovs{2},3)
	regTerm{1}(:,:,size(allcovs{1},3)+i) = allcovs{2}(:,:,i);
	regTerm{2}(:,:,size(allcovs{1},3)+i) = cspTerm{2};
end
opts.sym = 1;								%% use non-symmetric KL div for divCSP-WS
opts.lam = 0.1;							%% apply regularization
% opts.mode = 1;								%% maximize regularization term instead of minimizing
opts.beta = 0;							%% use beta divergence
% opts.pca = 0;
% opts.deflation = 1;	
% opts.beta = 1;
[P, results, y] = divcsp(cspTerm, regTerm, [], dim, opts);
P = P';


% covX1 = X1(:,:)-repmat(X1_mean,1,trial1);
% covX1 = covX1*covX1'/trial1;
% 
% covX2 = X2(:,:)-repmat(X2_mean,1,trial2);
% covX2 = covX2*covX2'/trial2;

% [V,D] = eig(covX1 ,covX1+ covX2);
% 
% r = abs(diag(D));
% [Y,I] = sort(r,'descend');
% V = V(:,I);
% r = Y;
% r = [r(1:n),r(chan1-n+1:chan1)];
% sf =[V(:,1:n),V(:,chan1-n+1:chan1)];
