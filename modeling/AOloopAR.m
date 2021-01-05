function [var_eps] = AOloopAR(G,H,A,phi,SNR,K)

n = size(H,1);      % dimension lifted wavefront  =(p+1)^2 
ns = size(G,1);     % dimension lifted sensor slopes  =2*p^2 
T = length(phi);   % number of temporal phase points  =(p+1)^2

% Initialiazation

epsk = zeros(n,T);  % residual wavefront
epskhat = zeros(n,T);
eps_piston_removed = zeros(n,T); % residual wavefront with mean removed
sk = zeros(ns,T);   % slopes measurements
u = zeros(n,T);     % input measurements 
sigma = zeros(T,1); % variance of the residual wavefront

% Matrices needed in the loop that can be computed once
M1  = (H'*H)\H'; M2 = A-K*G; M3 = A*H;

% Loop
for k = 1:(T-1)
    epsk(:,k+1) = phi(:,k+1)-H*u(:,k); % current residual wavefront
    sk(:,k+1) = awgn(G*epsk(:,k+1),SNR); % generation of the current 
                                         % slope measurement
    % Computation of the residual wavefront estimate needed to compute the 
    % appropriate input signal
    if k>1
        epskhat(:,k+1) = M2*epskhat(:,k)+K*sk(:,k)+M3*u(:,k-1)-H*u(:,k);
    else
        epskhat(:,k+1) = M2*epskhat(:,k)+K*sk(:,k)-H*u(:,k);
    end
    u(:,k+1) = M1*(M2*epskhat(:,k+1)+K*sk(:,k+1)+M3*u(:,k)); 
    % current input 
    eps_piston_removed(:,k+1) = epsk(:,k+1)-mean(epsk(:,k+1)); 
    % current residual wavefront without piston
    sigma(k+1) = var(eps_piston_removed(:,k+1)); % current variance
end

% variance of the residual wavefront
var_eps = mean(sigma);

end