

clc; close all; clear all;
global R22 Matrix XsN R32 S Phis1N V s  C A W R Phi0sN PhissN


load('systemMatrices.mat'); 
%matrices G(72x49)=(2p^2x(p+1)^2), 
%so p=6,H(49x49)=((p+1)^2x(p+1)^2),SNR=15
load('turbulenceData.mat'); %atmospheric turbulence
                            %data collected in open loop
                            %for 20 realizations
                            %matrix Phi (49x5000)
realizations = 20; % number of turbulence realizations
%% No Control
sigma_nocontrol = zeros(realizations,1);
for i = 1:realizations
    [sigma_nocontrol(i)] = AOloop_nocontrol(phiIdent{i},SNR,H,G);
end
figure
plot(1:realizations,sigma_nocontrol)
xlabel('Realizations');
ylabel('Variance of the residual wavefront')
title('No control action');

%% Static wavefront reconstruction
 sigma_e=0.01; %noise covariance
sigma_swr=zeros(realizations,1);
for i = 1:realizations
    C_phi0=(1/length(phiIdent{1}))*(phiIdent{i}*phiIdent{i}');
    [sigma_swr(i)] = AOloopMVM(G,H,C_phi0,sigma_e,phiIdent{i},SNR);
end
figure
plot(1:realizations,sigma_swr)
xlabel('Realizations');
ylabel('Variance of the residual wavefront')
title('Static Wavefront Reconstruction');

%% Vector Auto-Regressive model of Order 1
var_VAR=zeros(realizations,1);
for i = 1:realizations
    C_phi0=(1/length(phiIdent{1}))*(phiIdent{i}*phiIdent{i}');
    C_phi1=(1/(length(phiIdent{1})-1))*phiIdent{i}(:,2:end)*phiIdent{i}(:,1:(end-1))';
    [A,C_w,K] = computeKalmanAR(C_phi0,C_phi1,G,sigma_e);
    [var_VAR(i)] = AOloopAR(G,H,A,phiIdent{i},SNR,K);
end
figure
plot(1:realizations,var_VAR);
xlabel('Realizations');
ylabel('Variance of the residual wavefront')
title('Vector Auto-Regressive Model of Order 1');

% Subspace identification
Nid=3000;%number of points used for identification
Nval=2000;%number of points used for validation
n=[100 200 300 400 500 600 700 800 900 1000 1100]; %order of the system
s=30;

vaf=zeros(1,length(n));

phihat=phiIdent{1};
for i = 1:length(n)
    [~,~,~,vaf(i)]=n4sid(phihat,Nid,Nval,s,n(i));
end

figure()
plot(n,vaf,'xb', 'LineWidth', 2, 'MarkerSize', 12)
grid on
xlabel('model order')
ylabel('VAF')
title('VAF for different system orders and s=20')

%Matrices computation for the selected order
n_final=600;
s_final=20;

var_SI=zeros(realizations,1);
lamda=0.2;

for i = 1:realizations
phi=phiIdent{i};
[final_A,final_C,final_K,final_vaf] = n4sid(phi,Nid,Nval,s_final,n_final);
var_SI(i)=phiSid(H,final_A,final_K,final_C,lamda,phiSim{i});
end

figure()
plot(1:realizations,var_SI);
xlabel('Realizations');
ylabel('Variance of the residual wavefront')
title('Subspace identification using N4SID');