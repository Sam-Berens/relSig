%% RELSIG_EXAMPLE  Demonstration script for the relSig package
% Sam Berens (s.berens@sussex.ac.uk)
% RELSIG_EXAMPLE reproduces the workflow outlined in the README. It
% simulates data, fits an OLS model, performs family‑wise error‑rate
% control with RELSIG_FWE, and displays the significant results.
% 
% To run, simply type
%    relSig_Example
% at the MATLAB command prompt.
%
% See also OLS, RELSIG_FWE, BERENS_HOLM_SIDAK.
%
% -------------------------------------------------------------------------

%% Set a random seed
rng(0);

%% Set the ground thruth
ErrSigma = [...
    1.0,0.7,1.5;
    0.7,2.0,2.1;
    1.5,2.1,3.0];
n = 144;
X = [ones(n,1),randn(n,1)];
H = [0,1];
B = [...
    0.5,0.5,0.5;
    0.0,0.0,0.3];
Err = mvnrnd(zeros(1,3),ErrSigma,n);
Y = X*B + Err;

%% OLS estimation
[pValue,fStat,Bhat,C,Err] = OLS(Y,X,H);
ErrSigmaHat = cov(Err);

%% Run relSig
sig = relSig_FWE(ErrSigmaHat,X,H,pValue);