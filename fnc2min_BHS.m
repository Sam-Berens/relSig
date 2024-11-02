function [fVal] = fnc2min_BHS(mPrime,nullP)
Sig = false(size(nullP));
nIter = size(nullP,1);
for iIter = 1:nIter
    Sig(iIter,:) = berens_holm_sidak(nullP(iIter,:),mPrime);
end
fVal = (mean(any(Sig,2))-0.05).^2;
return