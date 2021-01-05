function [H]=compute_hankel(phi,i,s,N)
%start=0 for 1st matrix and start=s for future matrix
[n,Nt]=size(phi);%n=49, N=5000
identified_data=phi; %49*5000
identified_data=identified_data(:); %245000*1
H=zeros(s*n,Nt);
H(:,1)=identified_data((i*n)+1:(s*n+i*n));
for k=1:(Nt-s-i)
   h1=identified_data(((k+i)*n+1):(k+i+s)*n);
   H(:,k+1)=h1;
end
H=H(:,(1:N));
end