function [sig] = berens_holm_sidak(pValues,mPrime,alpha)

if nargin < 3
    alpha = 0.05;
end

%% Sort p-values in ascending order
[s_pValues,idx] = sort(pValues);
m = length(pValues); % Number of tests

%% Initialise the result
sig = false(size(pValues));

%% Berens-Holm-Šídák correction
kPrime = linspace(mPrime,1,m)';
for k = 1:m
    % Compute the Berens-Holm-Šídák adjusted alpha for the k-th test
    adjusted_alpha = 1 - (1 - alpha)^(1 / kPrime(k));
    
    % Check if the p-value is below the adjusted alpha threshold
    if s_pValues(k) <= adjusted_alpha
        % Mark as significant if the test passes
        sig(idx(k)) = true;
    else
        % Stop checking further if one test fails ...
        % ... (no further tests can be significant)
        break;
    end
end
return