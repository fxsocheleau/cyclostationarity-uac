%%
% This Matlab script  is an illustration of the MSML Doppler scale estimator described in Sec. III-C of the paper titled 
% "Cyclostationarity of Communication Signals in Underwater Acoustic Channels"
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
disp('================================================================================')
disp('-- Demo 1: CS-based Doppler scale estimator                                   --')
disp('-- QPSK signal in a MSML channel                                              --')         
disp('-- Channel simulated with Bellhop (2 paths with opposite-sign Doppler scales) --')
disp('================================================================================')

% 1- Load the QPSK signal (after MSML channel filtering)

load('MSML_Doppler_data')

% 2- Add noise

Eb_N0 = 5; % (Specify a given value for Eb/N0 in dB) 
Eb = sum(abs(MSML_sig_qpsk).^2)/nb_bits;
N0 = Eb*10.^(-Eb_N0/10);
noise = sqrt(N0/2)*(randn(size(MSML_sig_qpsk))+1i*randn(size(MSML_sig_qpsk)));
r = MSML_sig_qpsk + noise;

% 3- Apply the Doppler scale estimator (Alg. 1 of the paper)

v_max = 10; % maximum velocity in m/s
a_max = 1;  % maximum acceleration in m/s^2
c = 1500;   % assumed sound speed
[psi_hat,time,mu1_hat,mu2_hat] = f_estim_MSML_Doppler(r,fs,1/qpsk_rate,1+v_max/c,a_max/(2*c),D_max);


% 4- Display the results

disp(['Eb/N0 = ' num2str(Eb_N0) ' dB'])
disp(['Symbol rate = ' num2str(qpsk_rate) ' Bd'])
disp(['Duration = ' num2str(length(r)/fs) ' s'])

% NB: the estimation results may not be stored in the same order as the
% true values=> sort in ascending order for comparison.
% For a better interpretability, mu1_hat and mu2_hat are converted into
% velocity and acceleration with a sound speed arbitrarily set to 1500 m/s

[v_true_s,idxs_true]=sort(v_true); % sort the true velocities in ascending order
a_true_s = a_true(idxs_true);
v_hat = c*(1-mu1_hat); % sort the estimated velocities in ascending order
a_hat = -2*c*mu2_hat;
[v_hat_s,idxs_hat]=sort(v_hat);
a_hat_s = a_hat(idxs_hat);


for id = 1:D_max
    disp('---------')
    disp(['Scale #' num2str(id) ':'])
    disp('---------')
    disp(['Estimated velocity: ' num2str(v_hat_s(id),'%4.3f') ' m/s --> True velocity: ' num2str(v_true_s(id),'%4.3f') ' m/s'])
    disp(['Estimated acceleration: ' num2str(a_hat_s(id),'%4.3f') ' m/s^2 -->  True acceleration: ' num2str(a_true_s(id),'%4.3f') ' m/s^2'])
    disp(['Root-squared error (RSE) of the estimated time-warping function: ' num2str(10*log10(norm(psi_hat(:,idxs_hat(id))-psi_true(:,idxs_true(id)),2)./norm(psi_true(:,idxs_true(id)),2)),'%4.3f') ' dB'])
    
    % Plot the difference between the true and the estimated warping
    % functions
    figure
    plot(time, psi_true(:,idxs_true(id))-psi_hat(:,idxs_hat(id)),'Linewidth',2)
    grid on;
    hold on;
    plot(time, psi_true(:,id)-time,'Linewidth',2)
    legend('CSB','$\hat{\psi}(t)=t$','fontsize',18,'Interpreter','latex')
    xlabel('time (s)','Interpreter','latex','fontsize',18)
    ylabel('$\hat{\psi}(t)-\psi(t)$','Interpreter','latex','fontsize',18)
    title(['Time-warping function: estimation error, scale #' num2str(id)])
    
    saveas(gcf, ['../results/demo1_MSML_Doppler_fig' num2str(id) '.png'])
    
end
disp(' ');
disp('Figures are stored in the results folder')
