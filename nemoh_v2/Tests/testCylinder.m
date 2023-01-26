% This is a test function to test the workflow of NEMOH

clear
clc
close all




%% Step 1: Add required folders and make sure you are in the correct directory

% Add any required files
requiredfiles



% vhange current directory
[filepath,~,~] = fileparts(mfilename('fullpath'));
cd(filepath)



% Define some additional parameters in problem data structure (pdata)
pdata.plotdir      = mfoldername(mfilename('fullpath'),'plots');  % Plot directory
pdata.sname        = 'Cylinder_test';                             % Directory for storage of some results from Nemoh
pdata.meshfilename = 'axisym_new';


%% Step 2: Create a mesh 


n              = 3;                    % Number of points for discretization
radius         = 5;                    % Radius of the cylinder     
draft          = -10;                  % Draft  
r              = [radius radius 0];    % Radial coordinate in the polar system 
z              = [0 draft draft];      % Array of vertical coordinates



pdata.ntheta   = 50;                   % Number of angular discretization points
pdata.np       = 504;                  % Number of panels
pdata.zg       = 0;                    % Vertical position of gravity center
pdata.rho      = 1000;
pdata.g        = 9.81;
pdata.plotflag = 1;
pdata.meshflag = 1;
pdata.nbodies  = 1;


if pdata.meshflag
    [Mass,Inertia,KH,XB,YB,ZB] = axiMesh(r,z,n,0,0, 0,pdata);
end




%% Step 3: Run NEMOH

nbfreq         = 10;                                % number of calculations = # of BVP per freq * # of freq 
w              = linspace(0.04, 21,nbfreq)';         % periods from 3s to 20 s for waves
dir            = 0;                                  % Wave direction (angle of the incident wave)
depth          = 60;                                 % Water depth (0 for deep water)
tic
[A,B,Fe]       = Nemoh(w, dir, depth); 
toc


% Call Bemio functions for plotting and saving data into .h5 format
hydro         = struct();
dir_name      = mfoldername(mfilename('fullpath'),'Cylinder_test');
hydro         = readNEMOH(hydro,dir_name);
hydro         = radiationIRF(hydro,5,[],[],[],[]);
hydro         = radiationIRFSS(hydro,[],[]);
hydro         = excitationIRF(hydro,5,[],[],[],[]);
hydro.plotdir = pdata.plotdir;
writeBEMIOH5(hydro)

if pdata.plotflag
    plotBEMIO(hydro)
end
































%% Potential required files and folders

function requiredfiles
% Common files setup
% initialize
silentflag = 0; % silent

% add contents to path
RunSilent('AddSubmissionContents(mfilename)',silentflag)

% download required web zips
RunSilent('RequiredWebZips',silentflag)

% add contents to path (files have been downloaded)
RunSilent('AddSubmissionContents(mfilename)',silentflag)
end



function RequiredWebZips

disp('-> Obtaining required web zips')

% initialize index
ind = 0;

% initialize structure
zips = struct('url','','folder','','test','');

% zip
ind = ind + 1; % increment
zips(ind).url = 'https://github.com/AzadSaeed/Nemoh/archive/refs/heads/master.zip';
zips(ind).folder = 'NEMOH';
zips(ind).test = 'axiMesh';




% obtain full function path
full_fun_path = which(mfilename('fullpath'));
outputdir = fullfile(fileparts(full_fun_path),'lib');

% download and unzip
DownloadWebZips(zips,outputdir)

disp(' ')
end



function AddSubmissionContents(name)

disp('-> Adding submission contents to path')
disp(' ')

% turn off potential warning
warning('off','MATLAB:dispatcher:nameConflict')

% current file
fullfuncdir = which(name);

% current folder
submissiondir = fullfile(fileparts(fullfuncdir));

% add folders and subfolders to path
addpath(genpath(submissiondir))

% turn warning back on
warning('on','MATLAB:dispatcher:nameConflict')
end



function DownloadWebZips(zips,outputdir)

% store the current directory
olddir = pwd;

% create a folder for outputdir
if ~exist(outputdir, 'dir')
    mkdir(outputdir); % create the folder
else
    addpath(genpath(outputdir)); % add folders and subfolders to path
end

% change to the output directory
cd(outputdir)

% go through each zip
for k = 1:length(zips)

    % get data
    url = zips(k).url;
    folder = zips(k).folder;
    test = zips(k).test;

    % first check if the test file is in the path
    if exist(test,'file') == 0

        try
            % download zip file
            zipname = websave(folder,url);

            % save location
            outputdirname = fullfile(outputdir,folder);

            % create a folder utilizing name as the foldername name
            if ~exist(outputdirname, 'dir')
                mkdir(outputdirname);
            end

            % unzip the zip
            unzip(zipname,outputdirname);

            % delete the zip file
            delete([folder,'.zip'])

            % output to the command window
            disp(['Downloaded and unzipped ',folder])

        catch % failed to download
            % output to the command window
            disp(['Failed to download ',folder])

            % remove the html file
            delete([folder,'.html'])
        end

    else
        % output to the command window
        disp(['Already available ',folder])
    end
end

% change back to the original directory
cd(olddir)
end



function RunSilent(str,silentflag)

% if silent, capture the output
if silentflag
    O = evalc(str); %#ok<NASGU>
else
    eval(str);
end
end

