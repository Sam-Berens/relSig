function [sig] = relSig_FWE(ErrSigmaHat,X,H,pValue)

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