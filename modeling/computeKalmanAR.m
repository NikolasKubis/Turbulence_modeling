function [A,C_w,K] = computeKalmanAR(C_phi0,C_phi1,G,sigma_e)
A=C_phi1/C_phi0;
C_w=C_phi0-A*C_phi0*A';
R=(sigma_e^2)*eye(size(G,1));

% Check whether the conditions of Theorem 5.4 are met
if rank(obsv(A,G))<size(A,1)
    error('The (A,G) is not observable')
end
if rank(ctrb(A,C_w^(1/2)))<size(A,1)
    error('The (A,Cw^1/2) is not reachable')
end

[~,~,K] = dare(A',G',C_w,R);
K=K';
end