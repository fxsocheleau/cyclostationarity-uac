function J_DCS=f_signature_J_DCS(signal,fs,signature,fc,mu1_max,mu2_max,step_grid_mu1,step_grid_mu2)
%% J_DCS=f_signature_J_DCS(signal,fs,signature,fc,mu1_max,mu2_max,step_grid_mu1,step_grid_mu2)
%
% Compute the test statistic J_DCS (Eq. 56) in
% F.-X. Socheleau, "Cyclostationarity of Communication Signals in
% Underwater Acoustic Channels", IEEE JOE
%
% Input :  signal       : input signal to process
%          fs           : sampling frequency
%          signature    : matrix describing the cyclostationary signature (each row is
%          a pair (lag, cycle frequency)
%          fc           : carrier frequency
%          mu1_max      : maximum value for the first order coefficient in
%                         the optimization
%          mu2_max      : maximum value for the second order coefficient in
%                         the optimization
%          step_grid_mu1 : step of the grid search for mu1
%          step_grid_mu2 : step of the grid search for mu2
% Output:  
%          J_DCS      : test statistic
%
% Author: F.-X. Socheleau, IMT Atlantique, Lab-STICC, France
% Date: March 2022
%%


% Simple coarse-grid maximization

signal = signal(:);
t = (0:length(signal)-1).'/fs;
v_mu1 = [2-mu1_max:step_grid_mu1:mu1_max];
v_mu2 = [-mu2_max:step_grid_mu2:mu2_max];

J_DCS_grid = zeros(length(v_mu1),length(v_mu2));

for imu1=1:length(v_mu1) % loop over mu1
    mu1i = v_mu1(imu1);  
    for imu2=1:length(v_mu2) % loop over mu2
        mu2i = v_mu2(imu2);
        % dewarping
        epsilon_t = (mu1i-1)*t+mu2i*t.^2; 
        signalt=signal.*exp(-1i*2*pi*fc*epsilon_t);
        dwsignal = interp1(mu1i*t+mu2i*t.^2,signalt,t,'spline',0); 
        % compute CAF
        J_DCS_grid(imu1,imu2) = f_signature_J_CS(dwsignal,fs,signature);
    end
end
J_DCS = max(J_DCS_grid(:));
