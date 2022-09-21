function [psi_hat,time,mu1_hat,mu2_hat] = f_estim_MSML_Doppler(signal,fs,Ts,mu1_max,mu2_max,D_max)
%% [psi_hat,time,mu1_hat,m2_hat] = f_estim_MSML_Doppler(signal,fs,Ts,mu1_max,mu2_max,D_max)
%
% Estimate the time warping functions using Alg.1 described in the paper
% F.-X. Socheleau, "Cyclostationarity of Communication Signals in
% Underwater Acoustic Channels", IEEE JOE
%
% Input :  signal       : input signal to process
%          fs           : sampling frequency
%          Ts           : symbol period
%          mu1_max      : maximum value for the first order coefficient in
%                         the grid search
%          mu2_max      : maximum value for the second order coefficient in
%                         the grid search
%          D_max        : maximum number of Doppler scales
%        
% Output:  psi_hat      : matrix of size length(signal)xD_max
%          time         : time vector
%          mu1_hat      : estimated value of mu1 for each scale
%          mu2_hat      : estimated value of mu2 for each scale
%
% Author: F.-X. Socheleau, IMT Atlantique, Lab-STICC, France
% Date: March 2022
%%

%% Init outputs
psi_hat = zeros(length(signal),D_max);
time = (0:length(signal)-1).'/fs;
mu1_hat = zeros(D_max,1);
mu2_hat = zeros(D_max,1);
step_grid_mu1 = 1/9*1e-3; 
step_grid_mu2 = 1/6*1e-3;
%% Estimation
signal = signal(:);
r2 = abs(signal).^2;
for jd = 1:D_max % loop over the number of scales
    [mu1_hat(jd),mu2_hat(jd)]=f_maximize_cost_CS_Doppler(sqrt(r2),fs,1/Ts,mu1_max,mu2_max,step_grid_mu1,step_grid_mu2);
    psi_hat(:,jd)=mu1_hat(jd)*time+mu2_hat(jd)*time.^2;
    Gamma= mean(r2.*exp(-1i*2*pi*psi_hat(:,jd)/Ts));
    r2 = r2-2*real(Gamma.*exp(1i*2*pi*psi_hat(:,jd)/Ts));
end



