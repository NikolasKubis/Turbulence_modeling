function [var_eps] = AOloopMVM(G,H,C_phi0,sigma_e,phi,SNR)

%phi [(p+1)^2 x Nt]
%var_eps: variance of the residual wavefront
%Nt: time points

n = size(H,1);      % dimension lifted wavefront  =(p+1)^2 
ns = size(G,1);     % dimension lifted sensor slopes  =2*p^2
T = length(phi);   % number of temporal phase points  =(p+1)^2

epsk = zeros(n,T);  % residual wavefront
eps_piston_removed = zeros(n,T); % residual wavefront with mean removed
sk = zeros(ns,T);   % slopes measurements
u = zeros(n,T);     % input measurements 
sigma = zeros(T,1); % variance of the residual wavefront

% Matrices that we are gonna need for the wavefront reconstruction 
P = C_phi0;
W = (1/(sigma_e^2))*eye(ns);
M = ((H'*H)\H')*((inv(P)+G'*W*G)\G')*W;
%inv(H'*H)*H'*inv(inv(P)+G'*W*G)*G'*W;

for k = 1:(T-1)
    epsk(:,k+1) = phi(:,k+1)-H*u(:,k); % current residual wavefront
    eps_piston_removed(:,k+1) = epsk(:,k+1)-mean(epsk(:,k+1)); 
    sk(:,k+1) = awgn(G*epsk(:,k+1),SNR); % generation of 
                                         % the current slope measurement
    u(:,k+1) = M*sk(:,k+1)+u(:,k); % current input 
    sigma(k+1) = var(eps_piston_removed(:,k+1)); % current variance
end
var_eps = mean(sigma);

end