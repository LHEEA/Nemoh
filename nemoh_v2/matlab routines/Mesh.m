% 
% --> function [Mass,Inertia,KH,XB,YB,ZB]=Mesh(nBodies,n,X,tX,CG,nfobj)
%
% Purpose : Mesh generation of a symmetric body for use with Nemoh. Because
%           of symmetry, only half of the body must be described. 
%
% Inputs : description of body surface in large panels
%   - nBodies           : number of bodies
%   - n(nBodies)        : number of panels
%   - X(nBodies,n,4,3)  : coordinates of nodes of each panel
%   - tX(nBodies)       : translations
%   - CG(nBodies,3)     : position of gravity centre
%   - nfobj(nBodies)    : target number of panels for Aquaplus mesh
%
% Outputs : hydrostatics
%   - Mass(nBodies)         : masses of bodies
%   - Inertia(nBodies,6,6)  : inertia matrices (estimated assuming mass is distributed on
%   wetted surface)
%   - KH(nBodies,6,6)       : hydrostatic stiffness matrices
%   - XB,YB,ZB              : coordinates of buoyancy centers
%
% Copyright Ecole Centrale de Nantes 2014
% Licensed under the Apache License, Version 2.0
% Written by A. Babarit, LHEEA Lab.



function [Mass,Inertia,KH,XB,YB,ZB]=Mesh(nBodies,n,X,tX,CG,nfobj,data)

% Extract general information 
clear KH Mass Inertia XB YB ZB nx nf

rho     = data.rho;
g       = data.g;
nomrep  = data.Storage_dir;


% Create storage directory
system(['mkdir ',nomrep]);
system(['mkdir ',nomrep,filesep,'mesh']);
system(['mkdir ',nomrep,filesep,'results']);


% Initialize
Mass    = zeros(nBodies,1);
KH      = zeros(nBodies,6,6);
Inertia = zeros(nBodies,6,6);
XB      = zeros(nBodies,1);
YB      = zeros(nBodies,1);
ZB      = zeros(nBodies,1);
WPA     = zeros(nBodies,1);
nx      = zeros(nBodies,1);
nf      = zeros(nBodies,1);


% Symmetry check
sgn=X(1,1,1,2);
for c=1:nBodies
    for d=1:n(c)
        for i=1:4
            if (sgn*X(c,d,i,2) < 0)
                error('\n Be careful: it is assumed that a symmetry about the (xOz) plane is used. \n Only one half of the mesh must be descrided. \n The mesh will not be generated. \n');
                return;
            end
        end
    end
end



% Sauvegarde de la description du maillage
for c=1:nBodies
    fprintf('\n -> Meshing body number %g \n',c);
    clear x y z tri;
    fid=fopen([nomrep,filesep,'mesh',filesep,'mesh',int2str(c)],'w');
    fprintf(fid,'%g \n',4*n(c));
    fprintf(fid,'%g \n',n(c));
    nx(c)=0;
    for i=1:n(c)
        for j=1:4
            nx(c)=nx(c)+1;
            x(nx(c))=X(c,i,j,1);
            y(nx(c))=X(c,i,j,2);
            z(nx(c))=X(c,i,j,3);
            fprintf(fid,'%E %E %E \n',[X(c,i,j,1) X(c,i,j,2) X(c,i,j,3)]);
        end
    end
    for i=1:n(c)
        fprintf(fid,'%g %g %g %g \n',[4*(i-1)+1 4*(i-1)+2 4*(i-1)+3 4*(i-1)+4]');
    end
    fclose(fid);


    % Affichage de la description du maillage
    nftri=0;
    for i=1:n(c)
        nftri=nftri+1;
        tri(nftri,:)=[4*(i-1)+1 4*(i-1)+2 4*(i-1)+3];
        nftri=nftri+1;
        tri(nftri,:)=[4*(i-1)+1 4*(i-1)+3 4*(i-1)+4];
    end
    
    % Plot mesh panels
    if data.plotflag
        meshpanels(tri,x,y,z,nx,c,n,data);
    end

    fprintf('\n --> Number of nodes             : %g',nx(c));              % Added for consistency with v 3.01
    fprintf('\n --> Number of panels (max 2000) : %g \n',n(c));            % Added for consistency with v 3.01

    % Create Mesh.cal file
     CreateMeshcal(c,tX,CG,nfobj,rho,g);



%   Raffinement automatique du maillage et calculs hydrostatiques
    l = isunix;
    if l == 1
        system('mesh >Mesh.log');
    else
        %system(['mesh ',projdir]) % needs attention v 3.01
        filepath = which('Mesh.exe');
        system(filepath)
    end


%   Visualisation du maillage
    clear x y z NN nftri tri u v w xu yv zw;
    fid=fopen([nomrep,filesep,'mesh',filesep,'mesh',int2str(c),'.tec'],'r');
    fscanf(fid,'%s',2);
    nx(c)=fscanf(fid,'%g',1);
    fscanf(fid,'%s',2);
    nf(c)=fscanf(fid,'%g',1);
    fgetl(fid);
    fprintf('\n Characteristics of the mesh for Nemoh \n');
    fprintf('\n --> Number of nodes : %g',nx(c));
    fprintf('\n --> Number of panels : %g\n \n',nf(c));
    for i=1:nx(c)
        ligne=fscanf(fid,'%f',6);
        x(i)=ligne(1);
        y(i)=ligne(2);
        z(i)=ligne(3);
    end
    for i=1:nf(c)
        ligne=fscanf(fid,'%g',4);
        NN(1,i)=ligne(1);
        NN(2,i)=ligne(2);
        NN(3,i)=ligne(3);
        NN(4,i)=ligne(4);
    end
    nftri=0;
    for i=1:nf(c)
        nftri=nftri+1;
        tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
        nftri=nftri+1;
        tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
    end
    fgetl(fid);
    fgetl(fid);
    for i=1:nf(c)
        ligne=fscanf(fid,'%g %g',6);
        xu(i)=ligne(1);
        yv(i)=ligne(2);
        zw(i)=ligne(3);
        u(i)=ligne(4);
        v(i)=ligne(5);
        w(i)=ligne(6);
    end
    fclose(fid);

    if data.plotflag
        meshquiver(tri,x,y,z,xu,yv,zw,u,v,w,c,data);
    end


    % Read KH file
    KH = ReadKHfile(nomrep,c);


    % Read Hydrostatic data
    [Mass,Inertia]=ReadHydrodata(nomrep,c,rho,nBodies,data,CG);



%    To avoid some error when creating SM - Nan appears in place of CG
%     Cause of error is not knwon  
%     if data.B2B_SM
%         if c == 2
%             XB(c) = CG(c,1,1);
%             fscanf(fid,'%f',1)
%         else
%             XB(c)=fscanf(fid,'%f',1);
%         end
%     end

%    To make sure the saved file doesn't have the NaN
%     if data.B2B_SM
%         if  c ==2
%             fid=fopen([nomrep,filesep,'mesh',filesep,'Hydrostatics.dat'],'w');
%             fprintf(fid,'%g \n',CG(c,1,1));
%         end
%     end



%     if (~(c == nBodies))
%         next=input('Press enter to proceed with next body ');
%     end

%     if nBodies > 1
%         %fid=fopen([nomrep,filesep,'mesh',filesep,'Hydrostatics_',num2str(c-1),'.dat'],'w');
%         hydrofilename = strcat([nomrep,filesep,'mesh']);
%         filename = strcat('Hydrostatics_',num2str(c-1),'.dat');
%         cd(hydrofilename)
%         copyfile ("Hydrostatics.dat",filename)
%         filename = strcat('KH_',num2str(c-1),'.dat');
%         copyfile ("KH.dat",filename)
%         cd(data.mainpath)
%     end

end


% Write Nemoh input file
CreateNemohcalfile(nomrep,rho,g,nBodies,nx,nf,CG,c)


end





function meshpanels(tri,x,y,z,nx,c,n,data)
f = figure;
trimesh(tri,x,y,z,[zeros(nx(c),1)]);
title('Characteristics of the discretisation');
fprintf('\n --> Number of nodes             : %g',nx(c));
fprintf('\n --> Number of panels (max 2000) : %g \n',n(c));
xlabel('x')
ylabel('y')
zlabel('z')
figHandle = gca;
figname  = strcat('Mesh_Panels_',num2str(c));
filename = strcat(data.plot_dir,filesep,figname);
exportgraphics(figHandle,strcat(filename,'.pdf'));
close(f);
end





function CreateMeshcal(c,tX,CG,nfobj,rho,g)
fid=fopen('Mesh.cal','w');
fprintf(fid,['mesh',int2str(c),'\n'],1);
fprintf(fid,'1 \n %f 0. \n ',tX(c));
fprintf(fid,'%f %f %f \n',CG(c,:));
fprintf(fid,'%g \n 2 \n 0. \n 1.\n',nfobj(c));
fprintf(fid,'%f \n', rho);                        % Added for consistency with v 3.01
fprintf(fid,'%f \n', g);                          % Added for consistency with v 3.01
% fprintf(fid,'%f \n %f \n',[rho g]);             % Removed for consistency with v 3.01
status=fclose(fid);
end



function meshquiver(tri,x,y,z,xu,yv,zw,u,v,w,c,data)
    f = figure;
    trimesh(tri,x,y,z);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on;
    quiver3(xu,yv,zw,u,v,w);
    title('Mesh for Nemoh');
    figHandle = gca;
    figname  = strcat('Mesh_quiver_',num2str(c));
    filename = strcat(data.plot_dir,filesep,figname);
    exportgraphics(figHandle,strcat(filename,'.pdf'));
    close(f);
end



function KH = ReadKHfile(nomrep,c)
fid=fopen([nomrep,filesep,'mesh',filesep,'KH.dat'],'r');
for i=1:6
    ligne=fscanf(fid,'%g %g',6);
    KH(c,i,:)=ligne;
end
fclose(fid);
end






function [Mass,Inertia]=ReadHydrodata(nomrep,c,rho,nBodies,data,CG)
fid=fopen([nomrep,filesep,'mesh',filesep,'Hydrostatics.dat'],'r');
fscanf(fid,'%s',2);


%    To avoid some error when creating SM - Nan appears in place of CG
%     Cause of error is not knwon 
Err = 0;
try
    XB(c)=fscanf(fid,'%f',1);
catch ME
    warning('XB is NaN!');
    XB(c) = CG(c,1,1);
    Err = 1;
end

fgetl(fid);
fscanf(fid,'%s',2);
YB(c)=fscanf(fid,'%f',1);
fgetl(fid);
fscanf(fid,'%s',2);
ZB(c)=fscanf(fid,'%f',1);
fgetl(fid);
fscanf(fid,'%s',2);
Mass(c)=fscanf(fid,'%f',1)*rho;
fgetl(fid);
fscanf(fid,'%s',3);
WPA(c)=fscanf(fid,'%f',1);
fclose(fid);
clear ligne;
fid=fopen([nomrep,filesep,'mesh',filesep,'Inertia_hull.dat'],'r');
for i=1:3
    ligne=fscanf(fid,'%g %g',3);
    Inertia(c,i+3,4:6)=ligne;
end
Inertia(c,1,1)=Mass(c);
Inertia(c,2,2)=Mass(c);
Inertia(c,3,3)=Mass(c);
fclose(fid);

if Err == 1
    fid=fopen([nomrep,filesep,'mesh',filesep,'Hydrostatics.dat'],'r');
    fscanf(fid,'%s',23)
    Displacement = str2double(fscanf(fid,'%s',1));
    fscanf(fid,'%s',3)
    area = {};
    area{1} = fscanf(fid,'%s',1);
    
    T = readtable('Hydrostatics_0.dat','NumHeaderLines',0);
    T{1,3}= XB(c);
    T{1,7}= XB(c);
    T{4,3}= Displacement;
    T{5,4}= area;
    writetable(T,strcat(data.working_dir,filesep,data.Storage_dir,filesep,'mesh',filesep,'Hydrostatics.dat'))
end



%     fid=fopen([nomrep,filesep,'mesh',filesep,'Hydrostatics5.dat'],'w');
%     fscanf(fid,'%s',2);
%     fprintf(fid,'XF\n')
%     fprintf(fid,'=\n')
%     fprintf(fid,'%f',XB(c))
%     fclose(fid);





if nBodies > 1
    %fid=fopen([nomrep,filesep,'mesh',filesep,'Hydrostatics_',num2str(c-1),'.dat'],'w');
    hydrofilename = strcat([nomrep,filesep,'mesh']);
    filename = strcat('Hydrostatics_',num2str(c-1),'.dat');
    cd(hydrofilename)
    copyfile ("Hydrostatics.dat",filename)
    filename = strcat('KH_',num2str(c-1),'.dat');
    copyfile ("KH.dat",filename)
    cd(data.working_dir)
end





end






function CreateNemohcalfile(nomrep,rho,g,nBodies,nx,nf,CG,c)
fid=fopen([nomrep,filesep,'Nemoh.cal'],'w');
fprintf(fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
fprintf(fid,'%f				! RHO 			! KG/M**3 	! Fluid specific volume \n',rho);
fprintf(fid,'%f				! G			! M/S**2	! Gravity \n',g);
fprintf(fid,'0.                 ! DEPTH			! M		! Water depth\n');
fprintf(fid,'0.	0.              ! XEFF YEFF		! M		! Wave measurement point\n');
fprintf(fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n',c);
fprintf(fid,'%g				! Number of bodies\n',nBodies);
for c=1:nBodies
    fprintf(fid,'--- Body %g -----------------------------------------------------------------------------------------------------------------------\n',c);
    if isunix 
        fprintf(fid,['''',nomrep,filesep,'mesh',filesep,'mesh',int2str(c),'.dat''		! Name of mesh file\n']);
    else
        fprintf(fid,[nomrep,'\\mesh\\mesh',int2str(c),'.dat      ! Name of mesh file\n']);
    end
    fprintf(fid,'%g %g			! Number of points and number of panels 	\n',nx(c),nf(c));
    fprintf(fid,'6				! Number of degrees of freedom\n');
    fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Surge\n');
    fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Sway\n');
    fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Heave\n');
    fprintf(fid,'2 1. 0. 0. %f %f %f		! Roll about a point\n',CG(c,:));
    fprintf(fid,'2 0. 1. 0. %f %f %f		! Pitch about a point\n',CG(c,:));
    fprintf(fid,'2 0. 0. 1. %f %f %f		! Yaw about a point\n',CG(c,:));
    fprintf(fid,'6				! Number of resulting generalised forces\n');
    fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Force in x direction\n');
    fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Force in y direction\n');
    fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Force in z direction\n');
    fprintf(fid,'2 1. 0. 0. %f %f %f		! Moment force in x direction about a point\n',CG(c,:));
    fprintf(fid,'2 0. 1. 0. %f %f %f		! Moment force in y direction about a point\n',CG(c,:));
    fprintf(fid,'2 0. 0. 1. %f %f %f		! Moment force in z direction about a point\n',CG(c,:));
    fprintf(fid,'0				! Number of lines of additional information \n');
end
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


