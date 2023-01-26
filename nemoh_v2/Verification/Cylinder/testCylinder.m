% This is a test function to test the workflow of NEMOH

clear
clc
close all


%% Step 1: Create a mesh 

% Mesh generation for a vertical cylinder with 10 m diameter and 5 m draft
n  = 4;                   % Number of points for discretization
r  = [0, 5, 5, 0];        % Radial coordinate in the polar system 
z  = [5, 5, -5 -5];       % Array of vertical coordinates



data.ntheta   = 50;         % Number of angular discretization points
data.sname    = 'Cylinder'; % Directory name for storage of results
data.np       = 1000;         % Number of panels
data.zg       = -3;         % Vertical position of gravity center
data.rho      = 1025;
data.g        = 9.81;
data.plotflag = 1;

[Mass,Inertia,KH,XB,YB,ZB] = axiMesh(r,z,n, data);