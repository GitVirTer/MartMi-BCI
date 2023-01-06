function [ P, results, y ] = divcsp(X, Z, proj, dd, varargin)
%DIVCSP - Divegence-based Common Spatial Patterns
%
%Synopsis:
% [P, results, y] = divcsp(X, Z, proj, dd, opts)
%
%Arguments:
% X: cell containing average covariances of both classes, X{1} contains DxD cov mat of class 1, X{2} contains DxD cov mat of class 2
% Z: cell containing covariances for regularization, Z{1} and Z{2} contain DxDxn cov mats
% proj: ddxD projection matrix, only for testin purposes [] when computing spatial filters
% dd: dimensionality of extracted subspace
% varargin: optional parameters
%  max_iter: maximum number of iterations
%  nreps: number of repetitions
%  lam: lambda parameter in our regularization framework
%  beta: beta parameter
%  mode: if 0 substract regularization term, if 1 add regularization term
%  quiet: if 0 show output in each optimization step
%  deflation: if true use deflation algorithm, if false apply subspace method
%  csp_init: initialize with CSP rotation matrix
%  sym: use symmetric divergence in the regularization term
%  pca: apply PCA in the last step to obtain meaningful basis
%
%Output:
% P: Dxd projection matrix to the computed subspace
% results: data structure containing additional information
% y: objective function value
%
% Optimization algorithm adapted from
% von BÃ¼nau et al., Finding Stationary Subspaces in Multivariate Time Series, Phys Rev Let, 2009.
%
% Author(s): Wojciech Samek
% wojciech.samek@tu-berlin.de
%
% If using our code please cite one of the papers:
% W. Samek, M. Kawanabe, K.-R. MÃ¼ller, Divergence-based Framework for Common Spatial Patterns Algorithms
% IEEE Reviews in Biomedical Engineering, 2014.
%
% W. Samek, D. Blythe, K.-R. MÃ¼ller, M. Kawanabe, Robust Spatial Filtering with Beta Divergence
% Advances in Neural Information Processing Systems 26 (NIPS), 2013.
%

d = size(X{1}, 1);
opts = propertylist2struct(varargin{:});
opts = set_defaults(opts, ...
    	'max_iter', 100, ...
    	'nreps', 5, ...
    	'lam', 0, ...
    	'beta', 0, ...
    	'mode', 1, ...
		'quiet', false, ...
		'deflation', false, ...
		'csp_init', false, ...
		'sym', 1, ...
		'pca', false ...
    );

%% compute objective value and gradient for a fix projection, do not optimize
if(~isempty(proj))
	[y, grad] = objfun([], X, Z, dd, opts.mode, opts.lam, opts.beta, proj, opts.sym);
	P = proj;
	results.objfun_value = y;
	results.grad = grad;
	return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deflation mode 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.basis_proj = [];
if opts.deflation
	opts.deflation = false;
    opts.prev_y = []; %% Ìæ»»Î»ÖÃ
	opts_orig = opts;
	X_orig = X;
	w_init = [];
	
	while true 
		X = X_orig;
		opts = opts_orig;
		opts.init_rot = w_init;
		if isempty(opts.basis_proj) 
	  		[w, w_results, w_y] = divcsp(X, Z, proj, 1, opts);
			w_orth = w_results.rotation(2:end,:)*w_results.whitening;
			opts.basis_proj = eye(d);
		else
			opts.whitening = eye(d);
	  		[w, w_results, w_y] = divcsp(X, Z, proj, 1, opts);
			w_orth = w_results.rotation(2:end,:);
			w = w*opts.basis_proj;
		end
		if dd > 1
			for i=1:length(X)
				T{i} = zeros(d-1, d-1, size(X{i},3));
				for j=1:size(X{i},3)
					T{i}(:,:,j) = w_orth*X{i}(:,:,j)*w_orth';	
				end
			end
			X = T;
			for i=1:length(Z)
				U{i} = zeros(d-1, d-1, size(Z{i},3));
				for j=1:size(Z{i},3)
					U{i}(:,:,j) = w_orth*Z{i}(:,:,j)*w_orth';	
				end
			end
			Z = U;
			opts.prev_y = [ opts.prev_y w_results.objfun_value];
			opts.basis_proj = w_orth*opts.basis_proj;
			opts.deflation = true;
			[ P, P_results, P_y ] = divcsp(X, Z, proj, dd-1, opts);
		else
			P = []; P_y = []; P_results = {};
		end
			P = [ w; P ];
			y = [ w_y; P_y ];
			results = { w_results, P_results{:} };
			break;
	end
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over repetitions and return solution with lowest objective function value. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.nreps >= 1
  nreps = opts.nreps;
  r_P = cell(1, nreps);
  r_y = zeros(1, nreps);
  r_results = cell(1, nreps);
  opts.nreps = 0;
  for i=1:nreps
	 if(opts.csp_init==true && i==1)
	    [r_P{i}, r_results{i}, r_y(i) ] = divcsp(X, Z, proj, dd, opts);	% initialize with CSP solution
	 else
		 opts.csp_init = false;
	    [r_P{i}, r_results{i}, r_y(i) ] = divcsp(X, Z, proj, dd, opts);  % initialize with random rotation
	 end
  end
  [y, mini] = min(r_y);
  P = r_P{mini};
  results = r_results{mini};
  return; 
end

%% whitening matrix
W = inv(sqrtm(squeeze(mean(X{1},3) + mean(X{2},3))));
inv_W = sqrtm(squeeze(mean(X{1},3) + mean(X{2},3)));

%% rotation matrix
Binit = [];
if opts.csp_init == 1
	XX{1} = mult3(mean(X{1},3), W);
	[VV1 DD1] = eig(XX{1});
	[val ind] = sort([diag(DD1)' 1-diag(DD1)'],'descend');
	ind(ind > size(DD1,1)) = ind(ind > size(DD1,1)) - size(DD1,1);
	ind = ind(1:size(DD1,1));
	B = VV1(:,ind)';
else
	B = randrot(d);
end

ls_alpha = 0.5*(0.01+0.3);
ls_beta = 0.9;
direction = -1;
X_orig = X;
X{1} = mult3(X{1}, B*W);
X{2} = mult3(X{2}, B*W);
for ii=1:length(Z)
	Z{ii} = mult3(Z{ii}, B*W);
end
y_new = []; 
converged = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1;
for iter=1:opts.max_iter
   [y, grad, epoloss] = objfun([], X, Z, dd, opts.mode, opts.lam, opts.beta, proj, opts.sym);
	grad = grad*(-direction);
	y = y*(-direction);		
  if ~isempty(y_new) && y_new ~= y, error('Something is utterly wrong.\n'); end
  if ~opts.quiet, fprintf('iter=%d y=%.5g ||grad||=%.5g ', iter, y, norm(grad)); end
  if iter == 1
    alpha = -grad;
  else
    gamma = grad'*(grad-grad_old)/(grad_old'*grad_old);
    alpha = -grad;% + gamma*alpha_old;
  end
  grad_old = grad;
  alpha_old = alpha;
  alpha = alpha ./ (2*norm(alpha));
  M_alpha = reshape(alpha, [dd (d-dd)]);
  M_alpha = [ zeros(dd, dd) M_alpha; -M_alpha' zeros(d-dd, d-dd) ];
  for j=1:20
    M_new = t*M_alpha;
    R = expm(M_new);
 	 XX{1} = mult3(X{1}, R);
	 XX{2} = mult3(X{2}, R);
	 for ii=1:length(Z)
		ZZ{ii} = mult3(Z{ii}, R);
	 end
    y_new = objfun([], XX, ZZ, dd, opts.mode, opts.lam, opts.beta, proj, opts.sym);
	 y_new = y_new*(-direction);
    if y_new <= y %(y + ls_alpha*t*grad'*alpha)
      break;
    end
    t = ls_beta*t;
  end
  if y_new >= y
    if ~opts.quiet, fprintf('no step found\n'); end
    converged = true;
    break;
  end
  rel_dec = abs(y-y_new)/abs(y);
  rel_dec_thr = 1e-8;
  if rel_dec < rel_dec_thr
    if ~opts.quiet, fprintf('rel_dec < %f\n', rel_dec_thr); end
    converged = true;
    break;
  end
  if ~opts.quiet, fprintf('||step||=%.3g (%d) rel_dec=%.3g%%\n', t, j, 100*rel_dec); end
  X{1} = mult3(X{1}, R);
  X{2} = mult3(X{2}, R);
  for ii=1:length(Z)
	  Z{ii} = mult3(Z{ii}, R);
  end
  B = R*B;
end
if ~converged
  if ~opts.quiet, fprintf('Reached maximum number of iterations\n'); end
end
P = B(1:dd,:)*W;
if(dd > 1)
	if(opts.pca == true)
		[V D] = eig(P*mean(X_orig{1},3)*P');
		P = V'*P;	%% find meaningful basis in subspace by applying PCA
	end
end
results.objfun_value = y*(-direction);
results.whitening = W;
results.inv_whitening = inv_W;
results.rotation = B;
results.P = P;
results.iter = iter;
results.opts = opts;
results.A=inv_W*B';
results.epoloss = epoloss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = mult3(C, R)
[d1, d2, d3] = size(C);
C = reshape(C, [d1 d2*d3]);
C = reshape(R*C, [d1 d2 d3]);
C = permute(C, [2 1 3]);
C = reshape(C, [d1 d2*d3]);
C = reshape(R*C, [d1 d2 d3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fx, grad, epoloss] = objfun(M, C, Z, dd, mode, lam, beta, proj, sym)
[d,d,n] = size(C{1});
nn = size(Z{1},3);
if isempty(dd), dd = d; end
P = eye(dd,d);
if(~isempty(proj))
	P = proj;
end
if(~isempty(M))
	M = reshape(M, [dd (d-dd)]);
	M = [ zeros(dd, dd) M; -M' zeros(d-dd, d-dd) ];
	R = expm(M);
	P = R(1:dd,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CSP terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KLCSP term
if(beta == 0)
	for k=1:2
		pC{k} = zeros(dd,dd,n);
		inv_pC{k} = zeros(dd,dd,n);
		log_det{k} = zeros(1,n);
		trCC{k} = zeros(1, n);
		ipC{k} = zeros(dd, d, n);
		ipCC{k} = zeros(dd, d, n);
		ipCB{k} = zeros(dd, d, n);
		for i=1:n
		   pC{k}(:,:,i) = P*C{k}(:,:,i)*P';
			inv_pC{k}(:,:,i) = inv(pC{k}(:,:,i));
			log_det{k}(i) = log(det(pC{k}(:,:,i)));
		end
	end
	for i=1:n
	  trCC{1}(i) = trace(inv_pC{1}(:,:,i)*pC{2}(:,:,i));
	  trCC{2}(i) = trace(inv_pC{2}(:,:,i)*pC{1}(:,:,i));
	  ipC{1}(:,:,i) = inv_pC{1}(:,:,i)*P*C{1}(:,:,i);
	  ipC{2}(:,:,i) = inv_pC{2}(:,:,i)*P*C{2}(:,:,i);
	  ipCC{1}(:,:,i) = inv_pC{1}(:,:,i)*P*C{2}(:,:,i);
	  ipCC{2}(:,:,i) = inv_pC{2}(:,:,i)*P*C{1}(:,:,i);
	  ipCB{1}(:,:,i) = inv_pC{1}(:,:,i)*pC{2}(:,:,i)*inv_pC{1}(:,:,i)*P*C{1}(:,:,i);
	  ipCB{2}(:,:,i) = inv_pC{2}(:,:,i)*pC{1}(:,:,i)*inv_pC{2}(:,:,i)*P*C{2}(:,:,i);
   end
end

%% BetaCSP term
if(beta ~= 0)
	c1 = 1/((2*pi)^((beta*dd)/2) * (beta+1)^(dd/2));
	c2 = 1/((2*pi)^((beta*dd)/2));
	for k=1:2
		dett{k} = zeros(1, n);
		inv_pC{k} = zeros(dd, dd, n);
		pCBoth{k} = zeros(dd, dd, n);
		inv_pCBoth{k} = zeros(dd, dd, n);
		dettBoth{k} = zeros(1, n);
		for i=1:n
		  tmp = P*C{k}(:,:,i)*P';
		  dett{k}(i) = det(tmp);
		  inv_pC{k}(:,:,i) = inv(tmp);
		end
		ipC{k} = zeros(dd, d, n);
		ipCBoth{k} = zeros(dd, d, n);
	 end
	 CBoth{1} = beta*C{1} + C{2};
	 CBoth{2} = beta*C{2} + C{1};
	 for i=1:n
		pCBoth{1}(:,:,i) = P*CBoth{1}(:,:,i)*P';
		inv_pCBoth{1}(:,:,i) = inv(pCBoth{1}(:,:,i));
		dettBoth{1}(i) = det(pCBoth{1}(:,:,i));
		pCBoth{2}(:,:,i) = P*CBoth{2}(:,:,i)*P';
		inv_pCBoth{2}(:,:,i) = inv(pCBoth{2}(:,:,i));
		dettBoth{2}(i) = det(pCBoth{2}(:,:,i));
		ipC{1}(:,:,i) = inv_pC{1}(:,:,i)*P*C{1}(:,:,i);
		ipC{2}(:,:,i) = inv_pC{2}(:,:,i)*P*C{2}(:,:,i);
		ipCBoth{1}(:,:,i) = inv_pCBoth{1}(:,:,i)*P*CBoth{1}(:,:,i);
		ipCBoth{2}(:,:,i) = inv_pCBoth{2}(:,:,i)*P*CBoth{2}(:,:,i);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regularization terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(lam > 0)
	%% KL-nonstat
	if(beta == 0)
		for k=1:2
			zpC{k} = zeros(dd,dd,nn);
			zinv_pC{k} = zeros(dd,dd,nn);
			zlog_det{k} = zeros(1,nn);
			ztrCC{k} = zeros(1, nn);
			zipC{k} = zeros(dd, d, nn);
			zipCC{k} = zeros(dd, d, nn);
			zipCB{k} = zeros(dd, d, nn);
			for i=1:nn
				zpC{k}(:,:,i) = P*Z{k}(:,:,i)*P';
				zinv_pC{k}(:,:,i) = inv(zpC{k}(:,:,i));
				zlog_det{k}(i) = log(det(zpC{k}(:,:,i)));
			end
		end
		for i=1:nn
			ztrCC{1}(i) = trace(zinv_pC{1}(:,:,i)*zpC{2}(:,:,i));
			ztrCC{2}(i) = trace(zinv_pC{2}(:,:,i)*zpC{1}(:,:,i));
			zipC{1}(:,:,i) = zinv_pC{1}(:,:,i)*P*Z{1}(:,:,i);
			zipC{2}(:,:,i) = zinv_pC{2}(:,:,i)*P*Z{2}(:,:,i);
			zipCC{1}(:,:,i) = zinv_pC{1}(:,:,i)*P*Z{2}(:,:,i);
			zipCC{2}(:,:,i) = zinv_pC{2}(:,:,i)*P*Z{1}(:,:,i);
			zipCB{1}(:,:,i) = zinv_pC{1}(:,:,i)*zpC{2}(:,:,i)*zinv_pC{1}(:,:,i)*P*Z{1}(:,:,i);
			zipCB{2}(:,:,i) = zinv_pC{2}(:,:,i)*zpC{1}(:,:,i)*zinv_pC{2}(:,:,i)*P*Z{2}(:,:,i);
		end
	end

	%% Beta-nonstat
	if(beta ~= 0)
		for k=1:2
			zdett{k} = zeros(1, nn);
			zinv_pC{k} = zeros(dd, dd, nn);
			zpCBoth{k} = zeros(dd, dd, nn);
			zinv_pCBoth{k} = zeros(dd, dd, nn);
			zdettBoth{k} = zeros(1, nn);
			for i=1:nn
			  tmp = P*Z{k}(:,:,i)*P';
			  zdett{k}(i) = det(tmp);
			  zinv_pC{k}(:,:,i) = inv(tmp);
			end
			zipC{k} = zeros(dd, d, nn);
			zipCBoth{k} = zeros(dd, d, nn);
		 end
		 zCBoth{1} = beta*Z{1} + Z{2};
		 zCBoth{2} = beta*Z{2} + Z{1};
		 for i=1:nn
			zpCBoth{1}(:,:,i) = P*zCBoth{1}(:,:,i)*P';
			zinv_pCBoth{1}(:,:,i) = inv(zpCBoth{1}(:,:,i));
			zdettBoth{1}(i) = det(zpCBoth{1}(:,:,i));
			zpCBoth{2}(:,:,i) = P*zCBoth{2}(:,:,i)*P';
			zinv_pCBoth{2}(:,:,i) = inv(zpCBoth{2}(:,:,i));
			zdettBoth{2}(i) = det(zpCBoth{2}(:,:,i));
			zipC{1}(:,:,i) = zinv_pC{1}(:,:,i)*P*Z{1}(:,:,i);
			zipC{2}(:,:,i) = zinv_pC{2}(:,:,i)*P*Z{2}(:,:,i);
			zipCBoth{1}(:,:,i) = zinv_pCBoth{1}(:,:,i)*P*zCBoth{1}(:,:,i);
			zipCBoth{2}(:,:,i) = zinv_pCBoth{2}(:,:,i)*P*zCBoth{2}(:,:,i);
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpfx = 0; 
tmpgrad = zeros(dd,d);
if(mode == 0)
	sig = -1;
end
if(mode == 1)
	sig = 1;
end

%% KLCSP
if(beta == 0)
	epoloss.csp = 0.5*sum(-log_det{1} + log_det{2} + trCC{2} - dd)/n;
	epoloss.csp = epoloss.csp + 0.5*sum(-log_det{2} + log_det{1} + trCC{1} - dd)/n;
	tmpfx = tmpfx - (1-lam)*epoloss.csp;
	if(lam > 0)
		epoloss.reg = 0.5*sum(-zlog_det{1} + zlog_det{2} + ztrCC{2} - dd)/nn;
		if(sym == 1)
			epoloss.reg = epoloss.reg + 0.5*sum(-zlog_det{2} + zlog_det{1} + ztrCC{1} - dd)/nn;
		end
		tmpfx = tmpfx - sig*lam*epoloss.reg;
	end

   tmpgrad = tmpgrad - (1-lam)*0.5*sum(-2*ipC{1} + 2*ipC{2} - 2*ipCB{2} + 2*ipCC{2},3)/n;
	tmpgrad = tmpgrad - (1-lam)*0.5*sum(-2*ipC{2} + 2*ipC{1} - 2*ipCB{1} + 2*ipCC{1},3)/n;
	if(lam > 0)
		tmpgrad = tmpgrad - sig*lam*0.5*sum(-2*zipC{1} + 2*zipC{2} - 2*zipCB{2} + 2*zipCC{2},3)/nn;
		if(sym == 1)
			tmpgrad = tmpgrad - sig*lam*0.5*sum(-2*zipC{2} + 2*zipC{1} - 2*zipCB{1} + 2*zipCC{1},3)/nn;
		end
	end
end

%% BetaCSP
if(beta ~= 0)
	rdett{1} = reshape(repmat(dett{1},d*dd,1), [dd d n]);
	rdett{2} = reshape(repmat(dett{2},d*dd,1), [dd d n]);
	rdettBoth{1} = reshape(repmat(dettBoth{1},d*dd,1), [dd d n]);
	rdettBoth{2} = reshape(repmat(dettBoth{2},d*dd,1), [dd d n]);
	if(lam > 0)
		zrdett{1} = reshape(repmat(zdett{1},d*dd,1), [dd d nn]);
		zrdett{2} = reshape(repmat(zdett{2},d*dd,1), [dd d nn]);
		zrdettBoth{1} = reshape(repmat(zdettBoth{1},d*dd,1), [dd d nn]);
		zrdettBoth{2} = reshape(repmat(zdettBoth{2},d*dd,1), [dd d nn]);
	end

	epoloss.csp = sum((1/(beta*(beta+1))*c1*dett{1}.^(-beta/2) + 1/(beta+1)*c1*dett{2}.^(-beta/2) - 1/beta*c2*dett{2}.^(-(beta-1)/2).*dettBoth{1}.^(-1/2)))/n;
	epoloss.csp = 	epoloss.csp + sum((1/(beta*(beta+1))*c1*dett{2}.^(-beta/2) + 1/(beta+1)*c1*dett{1}.^(-beta/2) - 1/beta*c2*dett{1}.^(-(beta-1)/2).*dettBoth{2}.^(-1/2)))/n;
	tmpfx = tmpfx - (1-lam)*epoloss.csp;
	if(lam > 0)
		epoloss.reg = sum((1/(beta*(beta+1))*c1*zdett{1}.^(-beta/2) + 1/(beta+1)*c1*zdett{2}.^(-beta/2) - 1/beta*c2*zdett{2}.^(-(beta-1)/2).*zdettBoth{1}.^(-1/2)))/nn;
		if(sym == 1)
			epoloss.reg = epoloss.reg + sum((1/(beta*(beta+1))*c1*zdett{2}.^(-beta/2) + 1/(beta+1)*c1*zdett{1}.^(-beta/2) - 1/beta*c2*zdett{1}.^(-(beta-1)/2).*zdettBoth{2}.^(-1/2)))/nn;
		end
		tmpfx = tmpfx - sig*lam*epoloss.reg;
	end

	tmpgrad = tmpgrad - (1-lam)*sum(-1/(beta+1)*c1*rdett{1}.^(-beta/2).*ipC{1} - beta/(beta+1)*c1*rdett{2}.^(-beta/2).*ipC{2} + 1/beta*(beta-1)*c2*rdett{2}.^(-(beta-1)/2).*rdettBoth{1}.^(-1/2).*ipC{2} + 1/beta*c2*rdett{2}.^(-(beta-1)/2).*rdettBoth{1}.^(-1/2).*ipCBoth{1},3)/n;
	tmpgrad = tmpgrad - (1-lam)*sum(-1/(beta+1)*c1*rdett{2}.^(-beta/2).*ipC{2} - beta/(beta+1)*c1*rdett{1}.^(-beta/2).*ipC{1} + 1/beta*(beta-1)*c2*rdett{1}.^(-(beta-1)/2).*rdettBoth{2}.^(-1/2).*ipC{1} + 1/beta*c2*rdett{1}.^(-(beta-1)/2).*rdettBoth{2}.^(-1/2).*ipCBoth{2},3)/n;
	if(lam > 0)
		tmpgrad = tmpgrad - sig*lam*sum(-1/(beta+1)*c1*zrdett{1}.^(-beta/2).*zipC{1} - beta/(beta+1)*c1*zrdett{2}.^(-beta/2).*zipC{2} + 1/beta*(beta-1)*c2*zrdett{2}.^(-(beta-1)/2).*zrdettBoth{1}.^(-1/2).*zipC{2} + 1/beta*c2*zrdett{2}.^(-(beta-1)/2).*zrdettBoth{1}.^(-1/2).*zipCBoth{1},3)/nn;
		if(sym == 1)
			tmpgrad = tmpgrad - sig*lam*sum(-1/(beta+1)*c1*zrdett{2}.^(-beta/2).*zipC{2} - beta/(beta+1)*c1*zrdett{1}.^(-beta/2).*zipC{1} + 1/beta*(beta-1)*c2*zrdett{1}.^(-(beta-1)/2).*zrdettBoth{2}.^(-1/2).*zipC{1} + 1/beta*c2*zrdett{1}.^(-(beta-1)/2).*zrdettBoth{2}.^(-1/2).*zipCBoth{2},3)/nn;
		end
	end
end

grad = tmpgrad;
fx = tmpfx;
if nargout > 1
	grad = [ grad; zeros(d-dd, d) ];
	grad = grad - grad';
	grad = grad(1:dd,dd+1:end);
	grad = grad(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ R, M ] = randrot(d)
M = 10*(rand(d,d)-0.5);
M = 0.5*(M-M');
R = expm(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = propertylist2struct(varargin)
opt= [];
if nargin==0,
  return;
end

if isstruct(varargin{1}) | isempty(varargin{1}),
  opt= varargin{1};
  iListOffset= 1;
else
  iListOffset = 0;
end

nFields= (nargin-iListOffset)/2;
if nFields~=round(nFields),
  error('Invalid parameter/value list');
end

for ff= 1:nFields,
  fld = varargin{iListOffset+2*ff-1};
  if ~ischar(fld),
    error('Invalid parameter/value list');
  end
  opt.(fld)= varargin{iListOffset+2*ff};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt, isdefault]= set_defaults(opt, varargin)
isdefault= [];
if ~isempty(opt),
  for Fld=fieldnames(opt)',
    isdefault.(Fld{1})= 0;
  end
end

defopt = propertylist2struct(varargin{:});
for Fld= fieldnames(defopt)',
  fld= Fld{1};
  if ~isfield(opt, fld),
    [opt.(fld)]= deal(defopt.(fld));
    isdefault.(fld)= 1;
  end
end
