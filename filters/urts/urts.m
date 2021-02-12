function [M,P,D] = urts(M,P,a,dt,Q,param,alpha,beta,kappa,mat,same_p)
% urts_smooth1.m <- ORIGINAL
  if nargin < 5
    param = [];
  end
  if nargin < 6
    alpha = [];
  end
  if nargin < 7
    beta = [];
  end
  if nargin < 8
    kappa = [];
  end
  if nargin < 9
    mat = [];
  end
  if nargin < 10
    same_p = 1;
  end

  %
  % Apply defaults
  %
  if isempty(a)
    a = eye(size(M,1));
  end
  if isempty(Q)
    Q = zeros(size(M,1));
  end
  if isempty(mat)
    mat = 0;
  end

  %
  % Extend Q if NxN matrix
  %
  if size(Q,3)==1
    Q = repmat(Q,[1 1 size(M,2)]);
  end

  %
  % Run the smoother
  %
  D = zeros(size(M,1),size(M,1),size(M,2));
  for k=(size(M,2)-1):-1:1
    if isempty(param)
        params = [];
    elseif same_p
        params = param;
    else
        params = param{k};
    end
    [m_pred,P_pred,C] = ...
	ut_transform(M(:,k),P(:,:,k),a,dt,params,alpha,beta,kappa);
    P_pred = P_pred + Q(:,:,k);
    D(:,:,k) = C / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - m_pred);
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
  end