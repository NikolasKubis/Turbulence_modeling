function [phihat]=wavefront_reconstruction(G,C_phi0,sigma_e,phi,SNR)
ns = size(G,1);     % dimension lifted sensor slopes  =2*p^2
sk = awgn(G*phi,SNR);

phihat=C_phi0*G'*((G*C_phi0*G'+sigma_e^2*eye(ns))\sk);%phihat=C_phi0*G'*inv(G*C_phi0*G'+sigma_e^2*eye(ns))*sk;

end