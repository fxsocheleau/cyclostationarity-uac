function [C_hat,v_alpha, v_freq]  = f_spec_coherence_SIMO(signal,fs,v_alpha,nfft,noverlap,window)     
%% f_spec_coherence_SIMO(signal,fs,alpha_v,nfft,noverlap,window)     
%
% Compute the (cross-)spectral coherence density for SIMO signals (Eq. 71) in
% F.-X. Socheleau, "Cyclostationarity of Communication Signals in
% Underwater Acoustic Channels", IEEE JOE
%
% Input :  signal       : input signals to process
%          fs           : sampling frequency
%          v_alpha      : vector of cycle frequencies 
%          fc           : carrier frequency
%          nfft, noverlap, and window are as in function 'PWELCH' of Matlab.
% Output:  
%          C_hat        : estimated (cross-)coherence density of size
%          (size(signal,2) x size(signal,2) x nfft x length(v_alpha))
%          v_alpha      : vector of cycle frequencies 
%          v_freq       : vector of frequencies 
%
% Author: F.-X. Socheleau, IMT Atlantique, Lab-STICC, France
% Date: March 2022
%%

if length(window) == 1
    window = hanning(window);
end
window = window(:);

nwind = length(window);
% check inputs
v_alpha_n = v_alpha/fs;
if (max(v_alpha_n)>1)||(min(v_alpha_n)<0),error('alpha must be in [0,fs]'),end
if nwind <= noverlap,error('Window length must be > Noverlap');end
if nfft < nwind,error('Window length must be <= nfft');end

% init vectors and matrices
n_sig=size(signal,2);
C_hat = zeros(n_sig,n_sig,nfft,length(v_alpha_n));
n = size(signal,1);
t = (0:n-1)';
df = 1/nfft;
v_freq = (-1/2:df:1/2-df)*fs;
for isig1=1:n_sig
    y = signal(:,isig1);
    Sy = cyclic_spec(y,y,nfft,noverlap,window);
    for isig2=isig1:n_sig
        xorg = signal(:,isig2);
        if isig2==isig1 
            Sxref = Sy;
        else
            Sxref = cyclic_spec(xorg,xorg,nfft,noverlap,window);
        end
        for ialpha=1:length(v_alpha_n)
            alpha = v_alpha_n(ialpha);
            x = xorg.*exp(2i*pi*alpha*t);
            Syx = cyclic_spec(y,x,nfft,noverlap,window);
            Sx = circshift(Sxref,round(alpha/df)); % permutation = approximation to save time
            C_hat(isig1,isig2,:,ialpha) = Syx./sqrt(Sy.*Sx);
        end
    end
end
C_hat = fftshift(C_hat,3);
end


function [S,f]=cyclic_spec(x,y,nfft,noverlap,window)
% Estimate the cyclic spectrum using a periodogram
Lx = length(x);         % Number of data points
nwind = length(window); % length of window
K = fix((Lx-noverlap)/(nwind-noverlap));	% Number of windows
idx = 1:nwind;
f = (0:nfft-1)/nfft;
S = 0;
for i=1:K
    xw = window.*x(idx);
    yw = window.*y(idx);
    Yw = fft(yw,nfft);		
    Xw = fft(xw,nfft);		
    S = S+conj(Xw).*Yw;
    idx = idx + (nwind - noverlap);
end
end
 

