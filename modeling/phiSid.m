function [var_eps]=phiSid(H,A,K,C,lamda,phi)
[l,N] = size(phi);

[a,~]=size(A);

% Initialiazation
epsk = zeros(l,N);  % residual wavefront
eps_piston_removed = zeros(l,N); % residual wavefront with mean removed
u = zeros(l,N);     % input measurements 
sigma = zeros(N,1); % variance of the residual wavefront
x_sim=zeros(a,N);
phihat=zeros(l,N);% This is phihat: The estimation of the phi...

% Matrices needed in the loop that can be computed once
M1  = (H'*H+lamda*eye(size(H)))\H';
%M1  = inv(H'*H+lamda*eye(length(H)))*H'; 
 
 for k = 1:(N-1)
     
     x_sim(:,k+1)=A*x_sim(:,k)+K*(phi(:,k)-C*x_sim(:,k));
     % Model Simulation
     phihat(:,k+1)=C*x_sim(:,k+1);% phi prediction
     
     u(:,k) = M1*phihat(:,k+1);
     
     epsk(:,k+1) = phi(:,k+1)-H*u(:,k);
     eps_piston_removed(:,k+1) = epsk(:,k+1)-mean(epsk(:,k+1)); 
     % current residual wavefront without piston
     sigma(k+1) = var(eps_piston_removed(:,k+1)); % current variance 
 end
 % variance of the residual wavefront
var_eps = mean(sigma);
     
end