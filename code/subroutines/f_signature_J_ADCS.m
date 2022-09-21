function J_ADCS=f_signature_J_ADCS(signal,fs,signature,fc,mu1_max,mu2_max,step_grid_mu1,step_grid_mu2)
%% Synopsis: J_DCS=f_signature_J_ADCS(signal,fs,signature,fc,mu1_max,mu2_max,step_grid_mu1,step_grid_mu2)
%
% Compute the test statistic J_ADCS (Eq. 59) in
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
%          J_ADCS      : test statistic
%
% Author: F.-X. Socheleau, IMT Atlantique, Lab-STICC, France
% Date: March 2022
%%



% Init
signal = signal(:);
t = (0:length(signal)-1).'/fs;
v_mu1 = [2-mu1_max:step_grid_mu1:mu1_max];
v_mu2 = [-mu2_max:step_grid_mu2:mu2_max];
J_ADCS_grid = zeros(length(v_mu1),length(v_mu2));

% Simple coarse-grid maximization
for imu1=1:length(v_mu1)     % loop over mu1
    mu1i = v_mu1(imu1);  
    for imu2=1:length(v_mu2) % loop over mu2
        mu2i = v_mu2(imu2);
        % compute the test statistic for a given pair (mu1i,mu2i)
        idx_u_m1 = 0;
        for is=1:size(signature,1) % loop over signature
            u     = signature(is,1);
            alpha = signature(is,2);
            % compute the inverse warping function
            if mu2i~=0
                psi_inv_u =1/(2*mu2i)*(sqrt(mu1i^2+4*mu2i*u)-mu1i);
            else
                psi_inv_u = u/mu1i; 
            end
            idx_u = round(psi_inv_u*fs);
            if idx_u~=idx_u_m1
               corr_sig = conj(signal(1:end-idx_u)).*signal(idx_u+1:end); % correlation
               L_c = length(corr_sig);
            end
            corr_sig_s = corr_sig.*exp(-1i*2*pi*mu2i*2*u*fc*t(1:L_c));
            J_ADCS_grid(imu1,imu2) = J_ADCS_grid(imu1,imu2)+abs(mean(corr_sig_s.*exp(-1i*2*pi*(alpha*(mu1i*t(1:L_c)+mu2i*t(1:L_c).^2))))).^2; % cyclic correlation
            idx_u_m1 = idx_u;
        end

    end
end
J_ADCS = max(J_ADCS_grid(:));




