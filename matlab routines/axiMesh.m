% function [Mass,Inertia,KH,XB,YB,ZB]=axiMesh(r,z,n)
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

function [Mass,Inertia,KH,XB,YB,ZB]=axiMesh(r,z,n, data)


% Extract general information 
rho     = data.rho;
g       = data.g;
ntheta  = data.ntheta;
nomrep  = data.sname;
zG      = data.zg;
nfobj   = data.np;
theta   = 0.:pi/(ntheta-1):pi;
nx      = ntheta*n;
nf      = (n-1)*(ntheta-1);



% Calcul des sommets du maillage
x = kron(cos(theta),r);
y = kron(sin(theta),r);
z = repmat(z,1,ntheta);




% Calcul des facettes
I       = 1:n-1;
J       = 1:ntheta-1;
NN(1,:) = reshape((I(:)+n.*(J-1))',1,[]);
NN(2,:) = reshape((I(:)+1+n.*(J-1))',1,[]);
NN(3,:) = reshape((I(:)+1+n.*J)',1,[]);
NN(4,:) = reshape((I(:)+n.*J)',1,[]);



% Display mesh description
a_temp  = NN(1:3,:);
b_tempt = [NN(1,:);NN(3,:);NN(4,:)];
d_tempt = [a_temp;b_tempt];
d_tempt = reshape(d_tempt,3,[]);
tri     = d_tempt';



% Plot the triangular mesh data if plot flag is on
if data.plotflag
    figure;
    trimesh(tri,x,y,z,[zeros(nx,1)]);
    title('Characteristics of the discretisation');
    fprintf('\n --> Number of nodes             : %g',nx);
    fprintf('\n --> Number of panels (max 2000) : %g \n',nf);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    figHandle = gca;
    figname  = 'Mesh_Panels';
    filename = strcat(data.plotdir,figname);
    exportgraphics(figHandle,strcat(filename,'.pdf'));
end




% Create storage directory
system(['mkdir ',nomrep]);
system(['mkdir ',nomrep,filesep,'mesh']);
system(['mkdir ',nomrep,filesep,'results']);



% Create Mesh.cal file
fid = fopen('Mesh.cal','w');
fprintf(fid,'axisym \n',1);
fprintf(fid,'1 \n 0. 0. \n ');
fprintf(fid,'%f %f %f \n',[0. 0. zG]);
fprintf(fid,'%g \n 2 \n 0. \n 1.\n',nfobj);
fprintf(fid,'%f \n %f \n',[rho g]);
fclose(fid);



% Create and write to ID.dat file
fid=fopen([nomrep,filesep,'mesh',filesep,'axisym'],'w');
fprintf(fid,'%g \n',nx);
fprintf(fid,'%g \n',nf);
fprintf(fid,'%E %E %E \n', [x; y; z]);
fprintf(fid,'%g %g %g %g \n', NN);
fclose(fid);



% Run mesh.exe or mesh.log 
l = isunix;
if l == 1
    system('mesh >Mesh.log');
    % This part of if condition needs work
else
    filepath = which('Mesh.exe');
    system(filepath)
end



% Clear data
clear x y z NN nx nf a_temp b_temp d_temp tri u v w ;



% Read (get information about) axisym.tec file
fid   = fopen([nomrep,filesep,'mesh',filesep,'axisym.tec'],'r');
fscanf(fid,'%s',2);
nx    = fscanf(fid,'%g',1);
fscanf(fid,'%s',2);
nf    = fscanf(fid,'%g',1);
fgetl(fid);
fprintf('\n Characteristics of the mesh for Nemoh \n');
fprintf('\n --> Number of nodes : %g',nx);
fprintf('\n --> Number of panels : %g\n \n',nf);



% Read data from the open file into xyz_data matrix
xyz_data = fscanf(fid,'%f',[6,nx]);
x        = xyz_data(1,:);
y        = xyz_data(2,:);
z        = xyz_data(3,:);



% Read data from the open file into NN matrix
NN       = fscanf(fid,'%g',[4,nf]);



% Describe the triangular mesh
a_temp  = NN(1:3,:);
b_tempt = [NN(1,:);NN(3,:);NN(4,:)];
d_tempt = [a_temp;b_tempt];
d_tempt = reshape(d_tempt,3,[]);
tri     = d_tempt';
 


% Read additional data into vec_data matrix
fgetl(fid);
fgetl(fid);
vec_data   = fscanf(fid,'%g',[6,nf]);
xu         = vec_data(1,:);
yv         = vec_data(2,:);
zw         = vec_data(3,:);
u          = vec_data(4,:);
v          = vec_data(5,:);
w          = vec_data(6,:);
fclose(fid);


% Plot the mesh and the vector plot
if data.plotflag
    figure;
    trimesh(tri,x,y,z);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on;
    quiver3(xu,yv,zw,u,v,w);
    title('Mesh for Nemoh');
    figHandle = gca;
    figname  = 'Mesh_quiver';
    filename = strcat(data.plotdir,figname);
    exportgraphics(figHandle,strcat(filename,'.pdf'));
end


% Define the linear hydrodynamic restoring coefficient (or virtual mass)
clear KH;



% Open and read KH.dat file 
fid = fopen([nomrep,filesep,'mesh',filesep,'KH.dat'],'r');
KH = fscanf(fid,'%g %g',[6,6]);
fclose(fid);



% Clear variables
clear XB YB ZB Mass WPA Inertia a_temp b_tempt d_tempt



% Open and read Hydrostatics.dat file
Inertia=zeros(6,6);
fid=fopen([nomrep,filesep,'mesh',filesep,'Hydrostatics.dat'],'r');
fscanf(fid,'%s',2);
XB=fscanf(fid,'%f',1);
fgetl(fid);
fscanf(fid,'%s',2);
YB=fscanf(fid,'%f',1);
fgetl(fid);
fscanf(fid,'%s',2);
ZB=fscanf(fid,'%f',1);
fgetl(fid);
fscanf(fid,'%s',2);
Mass=fscanf(fid,'%f',1)*rho;
fgetl(fid);
fscanf(fid,'%s',2);
fscanf(fid,'%f',1);
fclose(fid);



% Open and read Inertia_hull.dat file
fid=fopen([nomrep,filesep,'mesh',filesep,'Inertia_hull.dat'],'r');
In       = fscanf(fid,'%g %g',[3,3]);
M        = Mass*eye(3);
Inertia  = blkdiag(M,In);
fclose(fid);



% Write Nemoh input file
fid=fopen([nomrep,filesep,'Nemoh.cal'],'w');
fprintf(fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
fprintf(fid,'%f				! RHO 			! KG/M**3 	! Fluid specific volume \n',rho);
fprintf(fid,'%f				! G			! M/S**2	! Gravity \n',g);
fprintf(fid,'0.                 ! DEPTH			! M		! Water depth\n');
fprintf(fid,'0.	0.              ! XEFF YEFF		! M		! Wave measurement point\n');
fprintf(fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n');
fprintf(fid,'1				! Number of bodies\n');
fprintf(fid,'--- Body 1 -----------------------------------------------------------------------------------------------------------------------\n');
if isunix 
    fprintf(fid,['''',nomrep,filesep,'mesh',filesep,'axisym.dat''		! Name of mesh file\n']);
else
    fprintf(fid,[nomrep,'\\mesh\\axisym.dat      ! Name of mesh file\n']);
end
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
fprintf(fid,'---');
fclose(fid);
fclose('all');

end

