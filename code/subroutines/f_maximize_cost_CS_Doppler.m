function [mu1_hat,mu2_hat]=f_maximize_cost_CS_Doppler(signal,fs,v_alpha,mu1_max,mu2_max,step_grid_mu1,step_grid_mu2)
%% [mu1_hat,mu2_hat]=f_maximize_cost_CS_Doppler(signal,fs,v_alpha,mu1_max,mu2_max,step_grid_mu1,step_grid_mu2)
%
% Maximize the cost function J (eq. 47) in
% F.-X. Socheleau, "Cyclostationarity of Communication Signals in
% Underwater Acoustic Channels", IEEE JOE
%
% Input :  signal       : input signal to process
%          fs           : sampling frequency
%          Ts           : symbol period
%          v_alpha      : values of cycle frequencies 
%          mu1_max      : maximum value for the first order coefficient in
%                         the grid search
%          mu2_max      : maximum value for the second order coefficient in
%                         the grid search
%          step_grid_mu1 : step of the grid search for mu1
%          step_grid_mu2 : step of the grid search for mu2
%        
% Output:  
%          mu1_hat      : estimated value of mu1 
%          mu2_hat      : estimated value of mu2
%
% Author: F.-X. Socheleau, IMT Atlantique, Lab-STICC, France
% Date: March 2022
%%


% Init search
signal = signal(:).';
t=[0:length(signal)-1]/fs;
v_mu1 = (2-mu1_max:step_grid_mu1:mu1_max);
v_mu2 = (-mu2_max:step_grid_mu2:mu2_max);

% Stopping criteria
C_EPSILON = 1e-10;
C_NITER = 30; % max number of iterations
stop=false;
iIter=1;
mu1i=0;
mu2i=0;
%% Maximization of the cost function
while iIter<C_NITER && stop==false
    if iIter==1  % first iteration: coarse grid search
        m_caf=compute_caf_grid_chirp_ZT(signal,v_mu1(1),step_grid_mu1, length(v_mu1),v_mu2,v_alpha,t);
        [Ji,ind]=max(m_caf(:));
        [idx_mu1,idx_mu2] = ind2sub(size(m_caf),ind);
        mu1i=v_mu1(idx_mu1);  
        mu2i=v_mu2(idx_mu2);
        mu1_hat=mu1i;
        mu2_hat=mu2i;
        Ji_max=Ji;
        v_grad_Ji=[0;0];
        mu1im1=mu1i-step_grid_mu1/16; % init guess
        mu2im1=mu2i-step_grid_mu2/8;
    else % next iterations: gradient ascent
        v_grad_Jim1=v_grad_Ji;
        v_grad_Ji=compute_caf_grad(signal,mu1i,mu2i,v_alpha,t);

        if norm(v_grad_Ji)/Ji<C_EPSILON
            stop=true;
        else
            mu1im2 = mu1im1;
            mu2im2 = mu2im1;
            mu1im1 = mu1i;
            mu2im1 = mu2i;
            kappai=abs(([mu1im1;mu2im1]-[mu1im2;mu2im2]).'*(v_grad_Ji-v_grad_Jim1))./(norm(v_grad_Ji-v_grad_Jim1)).^2; %Barzilai-Borwein step
            mu1i=mu1im1+kappai*v_grad_Ji(1);
            if mu2_max==0 
                mu2i = 0;
            else
                mu2i=mu2im1+kappai*v_grad_Ji(2);
            end

            Ji=compute_caf_grid(signal,mu1i,mu2i,v_alpha,t);
            if Ji > Ji_max % save results if better
               mu1_hat=mu1i; 
               mu2_hat=mu2i;
               Ji_max=Ji;
            end
        end
    end
    iIter=iIter+1;
end
%iIter
mu1_hat=min(mu1_hat,mu1_max);
mu1_hat=max(mu1_hat,-mu1_max);
mu2_hat=min(mu2_hat,mu2_max);
mu2_hat=max(mu2_hat,-mu2_max);
end

function m_caf=compute_caf_grid_chirp_ZT(signal,mu1_min,step_grid_mu1, Nmu1,v_mu2,v_alpha,t)
%% Compute the cost function using the chirp Z-transform
m_caf = zeros(Nmu1,length(v_mu2));
La=length(v_alpha);
Lx=length(signal);
fs=1/(t(2)-t(1)); % sampling frequency

for l=1:length(v_mu2)
    mu2=v_mu2(l);  
    for ia=1:La % loop over cycle frequencies 
        y=abs(signal).^2.*exp(-1i*2*pi*v_alpha(ia)*t.*(mu1_min+mu2*t));
        m_caf(:,l)=m_caf(:,l)+(abs(1/Lx*czt(y,Nmu1,exp(-1i*2*pi*v_alpha(ia)*(step_grid_mu1)/fs),1)).^2).';
    end
end
end

function m_caf=compute_caf_grid(signal,mu1,mu2,v_alpha,t)
%% Compute the cost function at a single point
m_caf = 0;
La=length(v_alpha);
for ia=1:La
    m_caf=m_caf+abs(mean(abs(signal).^2.*exp(-1i*2*pi*v_alpha(ia).*(mu1*t+mu2*t.^2)))).^2;             
end
end

function v_grad_Ji=compute_caf_grad(signal,mu1i,mu2i,v_alpha,t)
%% Compute the gradient of the cost function
psi_t = mu1i*t+mu2i*t.^2;
La=length(v_alpha);
v_grad_Ji=0;
for ia=1:La % loop over cycle frequencies 
        dv=mean(abs(signal).^2.*exp(-1i*2*pi*v_alpha(ia).*psi_t).*(-1i*2*pi*v_alpha(ia)*t));
        da=mean(abs(signal).^2.*exp(-1i*2*pi*v_alpha(ia).*psi_t).*(-1i*2*pi*v_alpha(ia)*t.^2));
        R=mean(abs(signal).^2.*exp(-1i*2*pi*v_alpha(ia).*psi_t));
        v_grad_Ji=v_grad_Ji+2*real(conj(R)*[dv;da]);
end

end