function [pValue,fStat,Bhat,C,Err] = OLS(Y,X,H)
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