function [A,C,K,vaf]=n4sid(phihat,Nid,Nval,s,n)
global R22 XsN R32 S Phis1N V W R Phi0sN PhissN

id_data = phihat(:,1:Nid); 
% data that will be used for identification
val_data = phihat(:,Nid+1:Nid+Nval); 
%data that will be used for validation

[l,N] = size(id_data);


Phi0sN = compute_hankel(id_data,0,s,(N-2*s+1));
PhissN = compute_hankel(id_data,s,s,(N-2*s+1));

% QR factorization
R = triu(qr([Phi0sN' PhissN']))';
R22 = R(1:s*l,1:s*l);
R32 = R((s*l)+1:2*s*l,1:s*l);


% page 332 book
[~,S,V]=svd((R32/R22)*Phi0sN); 
% lecture 6
%[~,S,V]=svd((1/N)*R32*R22');

SV=diag(S);
S = S(1:n,1:n);
V = V(:,1:n);
XsN = sqrtm(S)*V';

% figure()
% semilogy(SV, 'xb', 'LineWidth', 2, 'MarkerSize', 12);
% grid on
% xlabel('model order')
% ylabel('singular value')

columns=size(XsN);
Phis1N = compute_hankel(id_data,s,1,columns(2));

Y = [XsN(:,(2:end)) ; Phis1N(:,1:(end-1))];
F = XsN(:,1:(end-1));

AC = (Y*F')/(F*F');
A =AC(1:n,1:n);
C = AC(n+1:end,1:n);

W = XsN(:,(2:end)) - A*XsN(:,1:(end-1));
Vr = Phis1N(:,1:(end-1)) - C*XsN(:,1:(end-1));


Q = W*W';
Rr = Vr*Vr';
Sr = W*Vr';

[~, ~, K] = dare (A', C', Q, Rr, Sr);
K=K';

[a,~]=size(A);

[e,f]=size(val_data);

x_sim=zeros(a,f);
y_sim=zeros(e,f);

 % simulation of the system 
 for k=1:f
     x_sim(:,k+1)=A*x_sim(:,k)+K*(val_data(:,k)-C*x_sim(:,k));
     y_sim(:,k)=C*x_sim(:,k);
 end
 
 % VAF FINAL
num=0;
den=0;
for k=1:length(val_data)
   
    diff=(val_data(:,k)-y_sim(:,k));
    nr = norm(diff);
    nrsqr= nr^2;
    num = num + nrsqr;
    
    nr2 = norm(val_data(:,k));
    nrsqr2= nr2^2;
    den = den + nrsqr2;
    
end

num = num/length(val_data);
den = den/length(val_data);

vaf = max(0,(1-num/den)*100);


end