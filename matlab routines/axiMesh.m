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
%   - XB,YB,ZB  : coordintates of buoyancy center
%
% Warning : z(i) must be greater than z(i+1)
%
% Copyright Ecole Centrale de Nantes 2014
% Licensed under the Apache License, Version 2.0
% Written by A. Babarit, LHEEA Lab.
%
function [Mass,Inertia,KH,XB,YB,ZB]=axiMesh(r,z,n)
status=close('all');
ntheta=input('\n - Enter number of points for angular discretisation : ');
theta=[0.:pi/(ntheta-1):pi];
nx=0;
% Calcul des sommets du maillage
for j=1:ntheta
    for i=1:n    
        nx=nx+1;
        x(nx)=r(i)*cos(theta(j));
        y(nx)=r(i)*sin(theta(j));
        z(nx)=z(i);
    end;
end;
% Calcul des facettes
nf=0;
for i=1:n-1
    for j=1:ntheta-1
        nf=nf+1;
        NN(1,nf)=i+n*(j-1);
        NN(2,nf)=i+1+n*(j-1);
        NN(3,nf)=i+1+n*j;
        NN(4,nf)=i+n*j;
    end;
end;
% Affichage de la description du maillage
nftri=0;
for i=1:nf
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
end;
figure;
trimesh(tri,x,y,z,[zeros(nx,1)]);
title('Characteristics of the discretisation');
fprintf('\n --> Number of nodes             : %g',nx);
fprintf('\n --> Number of panels (max 2000) : %g \n',nf);
nomrep=input('\n - Directory name for storage of results : ');
system(['mkdir ',nomrep]);
system(['mkdir ',nomrep,'\Mesh']);
system(['mkdir ',nomrep,'\results']);
% Creation des fichiers de calcul du maillage
fid=fopen('mesh.cal','w');
fprintf(fid,'axisym \n',1);
fprintf(fid,'1 \n 0. 0. \n ');
zG=input(' - Vertical position of gravity center : ');
fprintf(fid,'%f %f %f \n',[0. 0. zG]);
nfobj=input(' - Target for number of panels : ');
fprintf(fid,'%g \n 2 \n 0. \n 1.\n',nfobj);
status=fclose(fid);
fid=fopen('ID.dat','w');
fprintf(fid,['% g \n',nomrep,' \n'],length(nomrep));
status=fclose(fid);
fid=fopen([nomrep,'/Mesh/axisym'],'w');
fprintf(fid,'%g \n',nx);
fprintf(fid,'%g \n',nf);
for i=1:nx
    fprintf(fid,'%E %E %E \n',[x(i) y(i) z(i)]);
end;
for i=1:nf
    fprintf(fid,'%g %g %g %g \n',NN(:,i)');
end;
status=fclose(fid);
% Raffinement automatique du maillage et calculs hydrostatiques
system('.\Mesh\Mesh.exe >Mesh\Mesh.log');
% Visualisation du maillage
clear x y z NN nx nf nftri tri u v w;
fid=fopen([nomrep,'\Mesh\axisym.tec'],'r');
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
end;
for i=1:nf
    ligne=fscanf(fid,'%g',4);
    NN(1,i)=ligne(1);
    NN(2,i)=ligne(2);
    NN(3,i)=ligne(3);
    NN(4,i)=ligne(4);
end;
nftri=0;
for i=1:nf
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
end;
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
end;
status=fclose(fid);
figure;
trimesh(tri,x,y,z);
hold on;
quiver3(xu,yv,zw,u,v,w);
title('Mesh for Nemoh');
clear KH;
KH=zeros(6,6);
fid=fopen([nomrep,'\Mesh\KH.dat'],'r');
for i=1:6   
    ligne=fscanf(fid,'%g %g',6);
    KH(i,:)=ligne;
end;
status=fclose(fid);
clear XB YB ZB Mass WPA Inertia
Inertia=zeros(6,6);
fid=fopen([nomrep,'\Mesh\Hydrostatics.dat'],'r');
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
Mass=fscanf(fid,'%f',1)*1025.;
ligne=fgetl(fid);
ligne=fscanf(fid,'%s',2);
WPA=fscanf(fid,'%f',1);
status=fclose(fid);
clear ligne
fid=fopen([nomrep,'\Mesh\Inertia_hull.dat'],'r');
for i=1:3
    ligne=fscanf(fid,'%g %g',3);
    Inertia(i+3,4:6)=ligne;
end;
Inertia(1,1)=Mass;
Inertia(2,2)=Mass;
Inertia(3,3)=Mass;
% Write Nemoh input file
fid=fopen([nomrep,'/Nemoh.cal'],'w');
fprintf(fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
fprintf(fid,'1000.0				! RHO 			! KG/M**3 	! Fluid specific volume \n');
fprintf(fid,'9.81				! G			! M/S**2	! Gravity \n');
fprintf(fid,'0.                 ! DEPTH			! M		! Water depth\n');
fprintf(fid,'0.	0.              ! XEFF YEFF		! M		! Wave measurement point\n');
fprintf(fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n');
fprintf(fid,'1				! Number of bodies\n');
fprintf(fid,'--- Body 1 -----------------------------------------------------------------------------------------------------------------------\n');
fprintf(fid,[nomrep,'\\mesh\\axisym.dat		! Name of mesh file\n']);
fprintf(fid,'%g %g			! Number of points and number of panels 	\n',nx,nf);
fprintf(fid,'6				! Number of degrees of freedom\n');
fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Surge\n');
fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Sway\n');
fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Heave\n');
fprintf(fid,'2 1. 0. 0. 0. 0. %f		! Roll about a point\n',zG);
fprintf(fid,'2 0. 1. 0. 0. 0. %f		! Pitch about a point\n',zG);
fprintf(fid,'2 0. 0. 1. 0. 0. %f		! Yaw about a point\n',zG);
fprintf(fid,'6				! Number of resulting generalised forces\n');
fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Force in x direction\n');
fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Force in y direction\n');
fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Force in z direction\n');
fprintf(fid,'2 1. 0. 0. 0. 0. %f		! Moment force in x direction about a point\n',zG);
fprintf(fid,'2 0. 1. 0. 0. 0. %f		! Moment force in y direction about a point\n',zG);
fprintf(fid,'2 0. 0. 1. 0. 0. %f		! Moment force in z direction about a point\n',zG);
fprintf(fid,'0				! Number of lines of additional information \n');
fprintf(fid,'--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'1	0.8	0.8		! Number of wave frequencies, Min, and Max (rad/s)\n');
fprintf(fid,'1	0.	0.		! Number of wave directions, Min and Max (degrees)\n');
fprintf(fid,'--- Post processing ---------------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'1	0.1	10.		! IRF 				! IRF calculation (0 for no calculation), time step and duration\n');
fprintf(fid,'0				! Show pressure\n');
fprintf(fid,'0	0.	180.		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n');
fprintf(fid,'0	50	400.	400.	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction\n');	
fprintf(fid,'---')
status=fclose(fid);
fclose('all');
end

