% 
% --> function [K, Mu]=FD2TD(w, A, B, T)
%
% Purpose:  Calculation of radiation coefficients in time domain from
%           frequency domain coefficients using Ogilvie's formulae
%
% Inputs :
% - w   : Frequency vector
% - A   : Added mass coefficient
% - B   : Damping coefficient
% - T   : Time vector
%
% Outputs :
% - K   : Retardation function K
% - Mu  : Added mass
%
% Copyright Ecole Centrale de Nantes 2014
% Licensed under the Apache License, Version 2.0
% Written by A. Babarit, LHEEA Lab.
%
function [K, Mu]=FD2TD(w, A, B, T)
clear K Mu;
n=length(T);
dt=T(2)-T(1);
nw=length(w);
dw=w(2)-w(1);
K=zeros(n,1);
for j=1:n
    for k=1:nw
    	K(j,1)=K(j,1)+B(k)*cos(w(k)*T(j))*dw;
    end;
end;
K=K*2./pi;
CM=zeros(nw,1);
Mu=0;
for k=1:nw
    for j=1:n
    	CM(k,1)=CM(k,1)+K(j)*sin(w(k)*T(j))*dt;
    end;
    CM(k,1)=A(k)+CM(k,1)/w(k);
    Mu=Mu+CM(k,1);
end;
Mu=Mu/nw;
for k=1:nw
    CM(k,1)=Mu;
end;
figure;
plot(w,B,'r',w,A,'b',w,CM,'b.');
xlabel('frequency (rad/s)');
ylabel('FD coefficents');
figure;
plot(T,K,'r');
xlabel('Time (s)');
ylabel('Retardation function');
%fid=fopen('Kt.dat','w');
%for j=1:n
%    fprintf(fid,'%E %E %E %E %E %E %E \n',T(j),0.,0.,0.,0.,K(j),0.);
%end;
%close(fid)
end