% 
% --> function [Mass,Inertia,KH,XB,YB,ZB]=axiMesh(r,z,n)
%
% Purpose : Mesh generation of an axisymmetric body for use with Nemoh
%
% Inputs : description of radial profile of the body
%   - n         : number of points for discretisation
%   - r         : array of radial coordinates
%   - z         : array of vertical coordinates
%
% Outputs : hydrostatics
%   - Mass      : mass of buoy
%   - Inertia   : inertia matrix (estimated assuming mass is distributed on
%   wetted surface)
%   - KH        : hydrostatic stiffness matrix
%   - XB,YB,ZB  : coordinates of buoyancy center
%
% Warning : z(i) must be greater than z(i+1)
%
% Copyright Ecole Centrale de Nantes 2014
% Licensed under the Apache License, Version 2.0
% Written by A. Babarit, LHEEA Lab.
%
function [Mass,Inertia,KH,XB,YB,ZB,nx,nf] = axiMesh_multi(r,z,n,in,num,coord,nomrep,mesh_file_name)     % SA

ind = num2str(num);
g = 9.81;
status=close('all');

theta=[0.:pi/(in.ntheta-1):pi];
nx=0;

% Calcul des sommets du maillage
for j=1:in.ntheta
    for i=1:n    
        nx=nx+1;
        x(nx)=r(i)*cos(theta(j)) + coord(1); % SA
        y(nx)=r(i)*sin(theta(j)) + coord(2); % SA
        z(nx)=z(i);
    end
end

% Calcul des facettes
nf=0;
for i=1:n-1
    for j=1:in.ntheta-1
        nf=nf+1;
        NN(1,nf)=i+n*(j-1);
        NN(2,nf)=i+1+n*(j-1);
        NN(3,nf)=i+1+n*j;
        NN(4,nf)=i+n*j;
    end
end

% Affichage de la description du maillage
nftri=0;
for i=1:nf
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
end

% figure;
% trimesh(tri,x,y,z,[zeros(nx,1)]);
% title('Characteristics of the discretisation');
% fprintf('\n --> Number of nodes             : %g',nx);
% fprintf('\n --> Number of panels (max 2000) : %g \n',nf);

% Creation des fichiers de calcul du maillage
fid=fopen('Mesh.cal','w');
fprintf(fid,mesh_file_name,1);
fprintf(fid,'\n 1 \n 0. 0. \n ');
fprintf(fid,'%f %f %f \n',[0. 0. in.zG]);
nfobj= 1; %input(' - Target for number of panels : ');
fprintf(fid,'%g \n 2 \n 0. \n 1.\n',nfobj);
fprintf(fid,'%f \n %f \n',[in.rho g]);

status=fclose(fid);
fid=fopen('ID.dat','w');
fprintf(fid,['% g \n',nomrep,' \n'],length(nomrep));

status=fclose(fid);
fid=fopen([nomrep,filesep,'mesh',filesep,mesh_file_name],'w');
fprintf(fid,'%g \n',nx);
fprintf(fid,'%g \n',nf);
for i=1:nx
    fprintf(fid,'%E %E %E \n',[x(i) y(i) z(i)]);
end

for i=1:nf
    fprintf(fid,'%g %g %g %g \n',NN(:,i)');
end

status=fclose(fid);
% Raffinement automatique du maillage et calculs hydrostatiques
l = isunix;
if l == 1
    system('mesh >Mesh.log');
else
    system('C:\Users\aksha\OneDrive\Documents\WEC_Project\WEC_Sim_NEMOH_Framework\WEC-Sim-4.4\Copy_of_WEC-Sim-4.4\new_model_test\Cylinder_test\hydroData\Mesh\Mesh.exe >C:\Users\aksha\OneDrive\Documents\WEC_Project\WEC_Sim_NEMOH_Framework\WEC-Sim-4.4\Copy_of_WEC-Sim-4.4\new_model_test\Cylinder_test\hydroData\Mesh\Mesh.log');
end

% Visualisation du maillage
clear x y z NN nx nf nftri tri u v w;
fid=fopen([nomrep,filesep,'mesh',filesep,horzcat(mesh_file_name,'.tec')],'r');
ligne=fscanf(fid,'%s',2);
nx=fscanf(fid,'%g',1);
ligne=fscanf(fid,'%s',2);
nf=fscanf(fid,'%g',1);
ligne=fgetl(fid);
fprintf('\n Characteristics of the mesh for Nemoh \n');
fprintf('\n --> Number of nodes : %g',nx);
fprintf('\n --> Number of panels : %g\n \n',nf);

for i=1:nx
    ligne=fscanf(fid,'%f',6);
    x(i)=ligne(1);
    y(i)=ligne(2);
    z(i)=ligne(3);
end

for i=1:nf
    ligne=fscanf(fid,'%g',4);
    NN(1,i)=ligne(1);
    NN(2,i)=ligne(2);
    NN(3,i)=ligne(3);
    NN(4,i)=ligne(4);
end

nftri=0;
for i=1:nf
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
end

ligne=fgetl(fid);
ligne=fgetl(fid);
for i=1:nf    
    ligne=fscanf(fid,'%g %g',6);
    xu(i)=ligne(1);
    yv(i)=ligne(2);
    zw(i)=ligne(3);
    u(i)=ligne(4);
    v(i)=ligne(5);
    w(i)=ligne(6);
end

status=fclose(fid);
figure;
trimesh(tri,x,y,z);
hold on;
quiver3(xu,yv,zw,u,v,w);
title('Mesh for Nemoh');

clear KH;
KH=zeros(6,6);
fid=fopen([nomrep,filesep,'mesh',filesep,'KH.dat'],'r');
for i=1:6   
    ligne=fscanf(fid,'%g %g',6);
    KH(i,:)=ligne;
end

status=fclose(fid);
clear XB YB ZB Mass WPA Inertia
Inertia=zeros(6,6);
fid=fopen([nomrep,filesep,'mesh',filesep,'Hydrostatics.dat'],'r');

ligne=fscanf(fid,'%s',2);
XB=fscanf(fid,'%f',1);
ligne=fgetl(fid);
ligne=fscanf(fid,'%s',2);
YB=fscanf(fid,'%f',1);
ligne=fgetl(fid);
ligne=fscanf(fid,'%s',2);
ZB=fscanf(fid,'%f',1);
ligne=fgetl(fid);
ligne=fscanf(fid,'%s',2);
Mass=fscanf(fid,'%f',1)*in.rho;
ligne=fgetl(fid);
ligne=fscanf(fid,'%s',2);
WPA=fscanf(fid,'%f',1);
status=fclose(fid);

clear ligne

fid=fopen([nomrep,filesep,'mesh',filesep,horzcat('Inertia_hull.dat')],'r');
for i=1:3
    ligne=fscanf(fid,'%g %g',3);
    Inertia(i+3,4:6)=ligne;
end

Inertia(1,1)=Mass;
Inertia(2,2)=Mass;
Inertia(3,3)=Mass;

fclose('all');

end
