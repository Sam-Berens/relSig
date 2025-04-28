function [sig] = relSig_FWE(ErrSigmaHat,X,H,pValue)
% RELSIG_FWE  Bootstrap‑calibrated FWE control for multivariate OLS tests.
% Sam Berens (s.berens@sussex.ac.uk)
% 
% SIG = RELSIG_FWE(ERRSIGMAHAT, X, H, PVALUE) returns a logical vector SIG
% indicating which raw p‑values PVALUE survive family‑wise error rate
% control at ALPHA = 0.05 using a bootstrap‑calibrated
% Berens‑Holm‑Šídák procedure.
% 
% The procedure:
%     1. Bootstraps the null distribution of p‑values by resampling from
%        a multivariate normal distribution with covariance ERRSIGMAHAT.
%     2. Estimates the effective number of independent tests (mPrime)
%        such that the empirical FWE equals ALPHA.
%     3. Applies BERENS_HOLM_SIDAK with that mPrime to the original
%        p‑values.
% 
% Syntax
% ------
% sig = relSig_FWE(ErrSigmaHat, X, H, pValue)
% 
% Input arguments
% ---------------
% ErrSigmaHat : m‑by‑m matrix
%               Estimated residual covariance matrix from OLS.
% X           : n‑by‑p design matrix used in the original model.
% H           : q‑by‑p contrast matrix defining the tested hypothesis
%               H * beta = 0.
% pValue      : 1‑by‑m vector of raw p‑values from OLS.
% 
% Output arguments
% ----------------
% sig : 1‑by‑m logical vector
%       TRUE for tests that survive FWE control.
% 
% Notes
% -----
% The number of bootstrap samples can be changed via the NITER constant
% near the top of this function.
% 
% Example
% -------
% See RELSIG_EXAMPLE for a complete demonstration.
% 
% See also OLS, BERENS_HOLM_SIDAK, FNC2MIN_BHS.
%
% -------------------------------------------------------------------------

%% Set the number of bootstrap samples
nIter = 1e5;

%% Bootstrap the null
m = numel(pValue);
n = size(X,1);
nullP = nan(nIter,m);
fh = waitbar(0,'Bootstrapping the null distribution...');
for iIter = 1:nIter
    nullY = mvnrnd(zeros(1,m),ErrSigmaHat,n);
    nullP(iIter,:) = OLS(nullY,X,H);
    if mod(iIter,97)==96
        waitbar(iIter/nIter,fh);
    end
end
close(fh);

%% Find mPrime
mPrime = fminbnd(...
    @(mm)fnc2min_BHS(mm,nullP),...
    1,m,optimset('Display','iter'));

%% Set the result
sig = berens_holm_sidak(pValue,mPrime);
return