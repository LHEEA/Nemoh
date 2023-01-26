% 
% --> function [A,B,Fe]=Nemoh(w, dir, depth)
%
% Purpose: Matlab wrapper for calculation of hydrodynamic coefficients using Nemoh
% 
% Inputs :
% - w     : Vector length(w) of wave frequencies (rad/s)
% - dir   : Wave direction (degrees)
% - depth : water depth (m), 0 for deep water.
%
% Outputs :
% - A  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of added mass coefficients
% - B  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of radiation damping coefficients
% - Fe : Matrix (6xnBodies)xlength(w) of exciation forces (complex
% values)
%
% Copyright Ecole Centrale de Nantes 2014
% Licensed under the Apache License, Version 2.0
% Written by A. Babarit, LHEEA Lab.

function [A,B,Fe]=Nemoh(w, dir, depth)

% Open ID.dat file
fid  = fopen('ID.dat');
fgetl(fid);
rep  = fscanf(fid,'%s');
fclose('all');



% Open and read Nemoh.cal file
fid=fopen([rep,filesep,'Nemoh.cal'],'r');
for i=1:6
    fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
fclose(fid);

% The following is slightly faster, but have to make sure it always works
% tic
% fid=fopen([rep,filesep,'Nemoh.cal'],'r');
% nBodies=textscan(fid,'%d',1,'delimiter','\n', 'headerlines',7);
% fclose(fid);
% t2 = toc



% Open and read Nemoh.cal file
fid=fopen([rep,filesep,'Nemoh.cal'],'r');
n=1;
clear textline;
textline={};
while (~feof(fid))
    textline(n)={fgetl(fid)};
    if (n == 4) 
        textline(n)={sprintf('%f                 ! DEPTH			! M		! Water depth',depth)};
    end
    if ((mod(n,18) == 9) && ((n-9)/18 <= nBodies))
        temp=cell2mat(textline(n));
        temp2=[];
        ntemp=length(temp);
        k=1;
        for i=1:ntemp
            if (temp(i) == '\')
                temp2=[temp2,temp(k:i),'\'];
                k=i+1;
            end            
        end
        temp2=[temp2,temp(k:ntemp)];
        textline(n)={temp2};
        cell2mat(textline(n));
    end
    if (n == 9+18*nBodies) 
        textline(n)={sprintf('%g %f %f       ! Number of wave frequencies, Min, and Max (rad/s)',length(w),w(1),w(length(w)))};
    end
     if (n == 10+18*nBodies) 
        textline(n)={sprintf('%g %f %f		! Number of wave directions, Min and Max (degrees)',1,dir,dir)};
    end
    n=n+1;
end
fclose(fid);


% Open and Write to Nemoh.cal file
fid = fopen([rep,filesep,'Nemoh.cal'], 'w'); 
for i=1:n-1
    fprintf(fid, [cell2mat(textline(i)),'\n']);
end
fclose(fid);


% Open and write to input.txt file in .txt format
fid=fopen([rep,filesep,'input.txt'],'wt');
fprintf(fid,' \n 0 \n');
fclose(fid);



% Calculate Hydrodynamic coefficients
l = isunix;
if l == 1
    fprintf('\n------ Starting NEMOH ----------- \n');
    system('preProc');
    fprintf('------ Solving BVPs ------------- \n');
    system('Solver');
    fprintf('------ Postprocessing results --- \n');
    system('postProc');
else
    fprintf('\n------ Starting NEMOH ----------- \n');
    %system('.\Nemoh\preProcessor.exe');
    filepath = which('preProcessor.exe');
    %filepath = 'A:\Azad\Resources\OpenWARP\source\NemohImproved\Nemoh\matlabRoutines\Nemoh\preProcessor.exe';
    system(filepath)
    fprintf('------ Solving BVPs ------------- \n');
    %system('.\Nemoh\Solver.exe');
    filepath = which('Solver.exe');
    %filepath = 'A:\Azad\Resources\OpenWARP\source\NemohImproved\Nemoh\matlabRoutines\Nemoh\Solver.exe';
    system(filepath)
    fprintf('------ Postprocessing results --- \n');
    %system('.\Nemoh\postProcessor.exe');
    filepath = which('postProcessor.exe');
    %filepath = 'A:\Azad\Resources\OpenWARP\source\NemohImproved\Nemoh\matlabRoutines\Nemoh\postProcessor.exe';
    system(filepath)
end
%% Lecture des resultats CA CM Fe
clear Periode A B Famp Fphi Fe;
fid=fopen([rep,filesep,'Nemoh.cal'],'r');
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
for i=1:2+18*nBodies
    ligne=fgetl(fid);
end
nw=fscanf(fid,'%g',1);
fclose(fid);
fid=fopen([rep,filesep,'results',filesep,'ExcitationForce.tec'],'r');
ligne=fgetl(fid);
for c=1:6*nBodies
    ligne=fgetl(fid);
end;
ligne=fgetl(fid);
for k=1:nw
    ligne=fscanf(fid,'%f',1+12*nBodies);
    w(i)=ligne(1);
    for j=1:6*nBodies
        Famp(k,j)=ligne(2*j);
        Fphi(k,j)=ligne(2*j+1);
    end;
end;
status=fclose(fid);
fid=fopen([rep,filesep,'results',filesep,'RadiationCoefficients.tec'],'r');
ligne=fgetl(fid);
for i=1:6*nBodies
    ligne=fgetl(fid);
end;
for i=1:nBodies*6
    ligne=fgetl(fid);
    for k=1:nw
        ligne=fscanf(fid,'%f',1+12*nBodies);
        for j=1:6*nBodies
            A(i,j,k)=ligne(2*j);
            B(i,j,k)=ligne(2*j+1);
        end;
        ligne=fgetl(fid);
    end;
end;
status=fclose(fid);
% Expression des efforts d excitation de houle sous la forme Famp*exp(i*Fphi)
i=sqrt(-1);
Fe=Famp.*(cos(Fphi)+i*sin(Fphi));
end