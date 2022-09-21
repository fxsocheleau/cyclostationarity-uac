%%
% This Matlab script  is an illustration of the cyclostationary-signature detectors described in the paper titled 
% "Cyclostationarity of Communication Signals in Underwater Acoustic
% Channels", IEEE JOE
% Author: F.-X. Socheleau
% IMT Atlantique, Lab-STICC, France. 
% March 2022
% Please make appropriate references to the corresponding paper if you use any of the matlab codes.
%%
clear all
close all


addpath(genpath(['.' filesep 'subroutines']))
addpath(genpath(['..' filesep 'data']))
disp(' ');
disp('=================================================')
disp('-- Demo 2: CS signature detection              --')
disp('-- CP-OFDM signal in a simulated MSML channel  --')  
disp('-- Impulsive noise (alpha-stable)              --')
disp('=================================================')

% 1- Load the OFDM signal (after MSML channel filtering) + signature data

load('MSML_CS_signature_data')

% 2- Compute the detection threshold for each detector
PFA = 1e-4;
threshold_J_CS      = f_threshold_pfa_cdf(PFA,x_cdf_CS,y_cdf_CS);
threshold_J_DCS     = f_threshold_pfa_cdf(PFA,x_cdf_DCS,y_cdf_DCS);
threshold_J_ADCS    = f_threshold_pfa_cdf(PFA,x_cdf_ADCS,y_cdf_ADCS);


% 2- Add impulsive noise
gSNR  = 5; % (Specify a value in dB for the generalized in-band SNR, NB: the variance of an alpha-stable noise does not exist for alpha<2 so we resort to another definition of SNR commonly used in the litterature ) 
Psig  = mean(abs(MSML_sig_ofdm).^2); % power of the signal
gamma = sqrt(Psig*fs/(2*B)*10^(-gSNR/10)); % scale parameter of the alpha-stable noise
alpha = 1.7;
noise = sqrt(1/2)*(random('Stable',alpha,0,gamma,0,size(MSML_sig_ofdm))+1i*random('Stable',alpha,0,gamma,0,size(MSML_sig_ofdm)));
r = MSML_sig_ofdm + noise;
% Plot the received signal
figure
plot((0:length(r)-1)/fs,real(r))
hold on;
grid on;
plot((0:length(r)-1)/fs,imag(r))
xlabel('time (s)')
legend('Re(r(t))', 'Im(r(t))')
title('Noisy OFDM signal')
saveas(gcf, ['../results/demo2_MSML_CS_signature_detection_fig1.png'])

% 3- Apply the three detectors described in Sec. III-D

s_r = f_spatial_sign(r); % apply the spatial-sign function on the observation
[J_CS]      =   f_signature_J_CS(s_r,fs,signature);                                                             % simple cyclostationarity-based test statistic
[J_DCS]     =   f_signature_J_DCS(s_r,fs,signature,fc,mu1_max,mu2_max,step_grid_mu1_DCS,step_grid_mu2_DCS);     % de-warped cyclostationarity-based test statistic
[J_ADCS]    =   f_signature_J_ADCS(s_r,fs,signature,fc,mu1_max,mu2_max,step_grid_mu1_ADCS,step_grid_mu2_ADCS);  % approximated de-warped cyclostationarity-based test statistic

% 4- Display the results
disp(['Pfa = ' num2str(PFA)])
disp(['Generalized SNR = ' num2str(gSNR) ' dB'])
disp(['Signature (lag,cycle frequency): ' mat2str(signature) ' (s,Hz)'])
disp('--------------------')
disp(' Detection results')
disp('--------------------')

if J_CS > threshold_J_CS
    disp('CS    => signature detected!')
else
    disp('CS    => signature NOT detected!')
end

if J_DCS > threshold_J_DCS
    disp('DCS   => signature detected!')
else
    disp('DCS   => signature NOT detected!')
end

if J_ADCS > threshold_J_ADCS
    disp('ADCS  => signature detected!')
else
    disp('ADCS  => signature NOT detected!')
end
disp(' ');
disp('Figures are stored in the results folder')
