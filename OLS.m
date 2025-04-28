function [pValue,fStat,Bhat,C,Err] = OLS(Y,X,H)
% OLS  Ordinary‑least‑squares estimation and F‑tests for linear contrasts.
% Sam Berens (s.berens@sussex.ac.uk)
% 
% [pValue, fStat, Bhat, C, Err] = OLS(Y, X, H) fits the multivariate
% linear model
%   Y = X * B + E
% where Y is n‑by‑m, X is n‑by‑p, and E is n‑by‑m random errors, and
% computes F‑statistics and p‑values for the null hypothesis
%   H * B = 0
% using the contrast matrix H (q‑by‑p).
% 
% Input arguments
% ---------------
% Y : n‑by‑m matrix
%     Response variables.
% X : n‑by‑p matrix
%     Design matrix (should include a column of ones for an intercept).
% H : q‑by‑p matrix
%     Contrast matrix defining the tested hypothesis.
% 
% Output arguments
% ----------------
% pValue : 1‑by‑m vector
%          Two‑sided p‑values for each response.
% fStat  : 1‑by‑m vector
%          Corresponding F‑statistics.
% Bhat   : p‑by‑m matrix
%          Least‑squares estimates of regression coefficients.
% C      : p‑by‑p‑by‑m array
%          Covariance of coefficient estimates for each response.
% Err    : n‑by‑m matrix
%          Residuals.
% 
% Example
% -------
% [p, ~, b] = OLS(Y, [ones(size(Y,1),1) X], H);
% 
% See also RELSIG_FWE.
% 
% -------------------------------------------------------------------------

Bhat = pinv(X)*Y;
P = X*pinv(X);
M = eye(size(P,1)) - P;
Err = M * Y;
dfMdl = rank(H);
dfErr = trace(M);
mse = diag((Err'*Err) / dfErr);
C = reshape(mse,1,1,size(Y,2)) .* pinv(X'*X);
fStat = nan(1,size(C,3));
for ii = 1:size(C,3)
    fStat(ii) = ...
        (H*Bhat(:,ii))' * (inv(H*C(:,:,ii)*H')) * (H*Bhat(:,ii)) / dfMdl;
end
pValue = 1-fcdf(fStat,dfMdl,dfErr);
return