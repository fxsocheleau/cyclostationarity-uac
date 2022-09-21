function [threshold,pfa_t] = f_threshold_pfa_cdf(pfa,x,cdf)
%% [threshold,pfa_t] = f_threshold_pfa_cdf(pfa,x,cdf)
%
% Compute the threshold to be used by a detector given a specified false
% alarm probability and the cumulative distribution function of the test
% statistic under H0
%
% Input :  pfa          : false-alarm probability (can be a vector)
%          x            : cumulative distribution function evaluated at x
%          cdf          : cumulative distribution function
%        
% Output:  
%          threshold  : threshold
%          pfa_t      : actual pfa estimated with the returned threshold 
% Author: F.-X. Socheleau, IMT Atlantique, Lab-STICC, France
% Date: March 2022
%%

if min(pfa)<min(1-cdf)
   error('PFA is too small')  
end
threshold=zeros(size(pfa));
for ia=1:length(pfa)
    % set a quantile (q) 
    q = 1-pfa(ia);
    pos = find(q<=cdf,1);
    threshold(ia) = x(pos);
    pfa_t=1-cdf(pos);
end
