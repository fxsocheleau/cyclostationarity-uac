function s_x = f_spatial_sign(x)
%% s_x = f_spatial_sign(x)
%
% Apply the spatial-sign function to an input signal x
%
% Input :  x   : input signal to process
%        
% Output:  
%          s_x : output signal
%
% Author: F.-X. Socheleau, IMT Atlantique, Lab-STICC, France
% Date: March 2022
%%
s_x         = zeros(size(x));
idx_nz      = find(x~=0);
s_x(idx_nz) = x(idx_nz)./abs(x(idx_nz));