function J_CS=f_signature_J_CS(signal,fs,signature) 
%% J_CS=f_signature_J_CS(signal,fs,signature) 
%
% Compute the test statistic J_CS (Eq. 55) in
% F.-X. Socheleau, "Cyclostationarity of Communication Signals in
% Underwater Acoustic Channels", IEEE JOE
%
% Input :  signal       : input signal to process
%          fs           : sampling frequency
%          signature    : matrix describing the cyclostationary signature (each row is
%          a pair (lag, cycle frequency)
%        
% Output:  
%          J_CS      : test statistic
%
% Author: F.-X. Socheleau, IMT Atlantique, Lab-STICC, France
% Date: March 2022
%%
signal = signal(:);
t = (0:length(signal)-1).'/fs;

J_CS = 0;
idx_u_m1 = 0;
for is=1:size(signature,1)
    u     = signature(is,1);
    alpha = signature(is,2);
    idx_u = round(u*fs);
    if idx_u~=idx_u_m1
        corr_sig = (conj(signal(1:end-idx_u)).*(signal(idx_u+1:end))); % correlation
        L_c = length(corr_sig);
    end
    J_CS = J_CS+abs(mean(corr_sig.*exp(-1i*2*pi*alpha*t(1:L_c)))).^2; % cyclic correlation
    idx_u_m1 = idx_u;
end


