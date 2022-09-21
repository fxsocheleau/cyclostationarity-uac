%%
% This Matlab script  is an illustration of the symbol-rate estimator described in the paper titled 
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
disp('====================================================')
disp('-- Demo 3: DSSS symbol-rate estimation            --')
disp('-- DSSS signal in a simulated dispersive channel  --')  
disp('====================================================')

% 1- Load the DSSS signal (after channel filtering) 

load('DISP_DSSS_data')

n_hydrophone = 1; % choose the number of hydrophones to use (must be less than 3)


% 2- Add noise

SNR     = 0; % (Specify a given value for the in-band SNR) 
Eb_N0   = SNR+10*log10(Ts*B/2); % deduce the value for Eb/N0;
signal  = DISP_DSSS_sig(:,1:n_hydrophone);
Eb      = sum(sum(abs(signal).^2,1),2)./(n_hydrophone*n_bits);
N0      = Eb*10.^(-Eb_N0/10);
noise   = sqrt(N0/2)*(randn(size(signal))+1i*randn(size(signal)));
r       = signal + noise;

switch n_hydrophone
    case 1
        disp('Using data from hydophone 1...')
        det_threshold =  det_threshold_h1;  
        figure
        [prr,f_dsp]=pwelch(r(:,1),1024,512,1024,fs,'centered');
        plot(f_dsp,10*log10(prr));
        xlabel('f (Hz)')
        ylabel('power/frequency (dB/Hz)')
        grid on;
        legend('Hydrophone #1')
        title('DSP of the received signal')
        xlim([-fs/2 fs/2])
    case 2    
        disp('Using data from hydophones 1 and 2...')
        det_threshold =  det_threshold_h2;
        figure
        [prr,f_dsp]=pwelch(r(:,1),1024,512,1024,fs,'centered');
        plot(f_dsp,10*log10(prr));
        grid on;
        hold on;
        [prr,f_dsp]=pwelch(r(:,2),1024,512,1024,fs,'centered');
        plot(f_dsp,10*log10(prr));
        xlabel('f (Hz)')
        ylabel('power/frequency (dB/Hz)')
        legend('Hydrophone #1','Hydrophone #2')
        title('DSP of the received signals')
        xlim([-fs/2 fs/2])
    otherwise
        disp('Using data from hydophones 1 to 3...')
        det_threshold =  det_threshold_h3;
        n_hydrophone = 3;
        figure
        [prr,f_dsp]=pwelch(r(:,1),1024,512,1024,fs,'centered');
        plot(f_dsp,10*log10(prr));
        grid on;
        hold on;
        [prr,f_dsp]=pwelch(r(:,2),1024,512,1024,fs,'centered');
        plot(f_dsp,10*log10(prr));
        [prr,f_dsp]=pwelch(r(:,3),1024,512,1024,fs,'centered');
        plot(f_dsp,10*log10(prr));
        xlabel('f (Hz)')
        ylabel('power/frequency (dB/Hz)')
        legend('Hydrophone #1','Hydrophone #2','Hydrophone #3')
        xlim([-fs/2 fs/2])
        title('DSP of the received signals')
        
end
saveas(gcf, ['../results/demo3_DISP_DSSS_symb_fig1.png'])

% 3- Compute the test statistic based on the coherence function

%  Parameters of the FAM estimator
NWind = 128;    % window size
nfft = 2*NWind; % size of fft
Noverlap = fix(2/3*NWind); %number of overlap samples 
alphamin = 0; % smallest cycle frequency
alphamax = B; % largest cycle frequency
L = size(DISP_DSSS_sig,1);		% signal length
da = fs/L;           % cyclic frequency resolution
v_alpha = (alphamin:da:alphamax);
% compute the spectral coherence density
disp('Computation of the spectral coherence density...')
C_hat = f_spec_coherence_SIMO(r,fs,v_alpha,nfft,Noverlap,NWind);



% 3- Detect the significant cycle frequencies

test_stat = max(abs(squeeze(sum(sum(C_hat,1),2))));

figure % plot test statistic
plot(v_alpha,test_stat)
xlim([0 alphamax])
grid on 
hold on
xlabel('\alpha (Hz)')
ylabel('test statistic')


test_stat(test_stat<det_threshold)=0;
[pks,Aset_hat_idx] = findpeaks(test_stat); % estimated cycle frequencies
Aset_hat = (Aset_hat_idx-1)*da;
plot(Aset_hat,pks,'r+') % plot the detected cycle frequencies
title('Detected cycle frequencies')
saveas(gcf, ['../results/demo3_DISP_DSSS_symb_fig2.png'])
% Possible improvement (compared to the article): set amplitude-dependent weights to each candidate cycle frequency 
norm_set = max(pks)/100;
weights = round(pks(:)/norm_set);
Aset_hat_w = ones(1,sum(weights));
idxc = 1;
for ipks = 1:length(pks)
    Aset_hat_w(idxc:idxc+weights(ipks)-1)=Aset_hat(ipks);
    idxc = idxc+weights(ipks);
end

% 4- Find the most likely  greatest common divisor (GCD) of the estimated cycle frequencies
Nc_max = 64; % rough bounds on the number of chips 
Nc_min = 2;
epsilon = 0.05;
if ~isempty(Aset_hat_w)
    Aset_hat_w = [0 Aset_hat_w];
    dAset_hat = abs(bsxfun(@minus,Aset_hat_w,Aset_hat_w'));
    Dset_hat = unique(dAset_hat(:));
    Dset_hat(Dset_hat<B/Nc_max | Dset_hat>B/Nc_min)=[];
    count_max = 0;
    numel = zeros(length(Dset_hat),1);
    for ida=1:length(Dset_hat)
        da_test = Dset_hat(ida);
        numel(ida)=length(find(mod(Aset_hat_w,da_test)<=da_test*epsilon | mod(Aset_hat_w,da_test)>=da_test*(1-epsilon)));
    end

    figure % plot counting results
    p=plot(Dset_hat,numel,'+-');
    p.Color = [0.8500 0.3250 0.0980];
    grid on;
    xlabel('\Delta_{\alpha(d)} (Hz)')
    ylabel('c(d)')
    title('Counting results')
    saveas(gcf, ['../results/demo3_DISP_DSSS_symb_fig3.png'])
    [~,idxm]=max(numel);
    Ts_hat = 1./Dset_hat(idxm);
else
    Ts_hat = Inf;   
end

% 4- Display the results
disp(['SNR = ' num2str(SNR) ' dB'])
disp(['Carrier frequency = ' num2str(fc) ' (Hz)'])
disp(['Bandwidth = ' num2str(B) ' (Hz)'])
disp(['Duration = ' num2str(size(r,1)/fs) ' (s)'])
disp('--------------------')
disp(' Estimation results')
disp('--------------------')
disp(['Estimated symbol-rate: ' num2str(1/Ts_hat,'%4.2f') ' (Hz)'])
disp(['=> True symbol-rate: ' num2str(1/Ts,'%4.2f') ' (Hz)'])
disp(['Estimated symbol duration: ' num2str(Ts_hat,'%4.3f') ' (s)'])
disp(['=> True symbol duration: ' num2str(Ts,'%4.3f') ' (s)'])
disp(['Relative estimation error: ' num2str(abs(Ts_hat-Ts)/Ts*100,'%4.2f') ' %'])

if SNR>= 10
    disp('--------------------------------------')
    disp('Warning: for high SNRs, it may be necessary to increase the detection threshold and/or to limit the computation of the spectral coherence density to the support [-B+\alpha;B], where B denotes the monolateral bandwidth.')
end

disp(' ');
disp('Figures are stored in the results folder')


