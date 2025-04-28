function [sig] = berens_holm_sidak(pValues,mPrime,alpha)
% Berens-Holm-Sidak step‑down multiple‑testing correction.
% Sam Berens (s.berens@sussex.ac.uk)
% 
% SIG = BERENS_HOLM_SIDAK(PVALUES, MPRIME) returns a logical vector SIG
% whose TRUE entries mark hypotheses that remain significant after the
% Berens‑Holm‑Šídák step‑down procedure controlling the family‑wise error
% rate (FWE) at ALPHA = 0.05.
% 
% SIG = BERENS_HOLM_SIDAK(PVALUES, MPRIME, ALPHA) uses the specified FWE
% threshold ALPHA (0 < ALPHA < 1).
% 
% Input arguments
% ---------------
% pValues : vector
%           Raw p‑values (one per hypothesis).
% mPrime  : scalar
%           Effective number of independent tests.  See RELSIG_FWE for how 
%           this is estimated.
% alpha   : scalar, optional (default 0.05).
%           Target family‑wise error rate.
%
% Output arguments
% ----------------
% sig     : logical vector (same size as pValues)
%           TRUE for hypotheses that remain significant after correction.
% 
% Example
% -------
% sig = berens_holm_sidak([0.01 0.20 0.03], 2);
% 
% See also RELSIG_FWE, FNC2MIN_BHS.
% 
% -------------------------------------------------------------------------

%% Check inputs
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