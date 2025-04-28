function [fVal] = fnc2min_BHS(mPrime,nullP)
% FNC2MIN_BHS  Objective function for calibrating Berens‑Holm‑Šídák FWE.
% Sam Berens (s.berens@sussex.ac.uk)
% 
% fVal = fnc2min_BHS(mPrime, nullP) returns the squared error between the
% empirical family‑wise error rate obtained from applying
% BERENS_HOLM_SIDAK with parameter mPrime to the bootstrap null p‑values
% in NULLP and the nominal rate ALPHA = 0.05.
% 
% The function is intended for use with optimisation routines such as
% FMINBND to find the value of mPrime that yields the desired FWE.
% 
% Input arguments
% ---------------
% mPrime : scalar
%          Candidate effective number of independent tests.
% nullP  : nIter‑by‑m matrix
%          Bootstrap p‑values obtained under the null hypothesis.
% 
% Output arguments
% ----------------
% fVal   : scalar
%          Squared difference between empirical and nominal FWE.
% 
% See also BERENS_HOLM_SIDAK, RELSIG_FWE.
% 
% -------------------------------------------------------------------------

Sig = false(size(nullP));
nIter = size(nullP,1);
for iIter = 1:nIter
    Sig(iIter,:) = berens_holm_sidak(nullP(iIter,:),mPrime);
end
fVal = (mean(any(Sig,2))-0.05).^2;
return