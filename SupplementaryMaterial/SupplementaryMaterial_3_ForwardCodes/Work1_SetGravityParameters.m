% This Matlab script shows how to generate the parameter files "TFG.ForPar"
% needed for the three cpp main functions:
% Work1_GravityEstimate.cpp, Work1_GravityEstimateOpenMP.cpp, and Work1_GravityEstimateMPI.cpp.
% For more details, see Zhang P. et al.,
% "An alternative approach for accurately calculating gravitational and magnetic fields of a spherical prism."

% For the components of the binary parameter file, please refer to line 44

%%

clc;
clear;
close all;

%%

Fai1=[30;30;40;40];
Fai2=[40;40;50;50];
Lamda1=[110;120;110;120];
Lamda2=[120;130;120;130];

NumDensityModel=size(Fai1,1);

RE=6371.2;

R1=(RE-30).*ones(NumDensityModel,1);
R2=RE.*ones(NumDensityModel,1);
rho=1.*ones(NumDensityModel,1);

%%

AbsTol=1e-16;
RelTol=1e-16;

%%

[Longitude,Latitude]=meshgrid(100.05:0.1:140.05,20.05:0.1:60.05);
Longitude=Longitude(:);
Latitude=Latitude(:);

NumCoordinates=size(Longitude,1);
Radius=(RE+1).*ones(NumCoordinates,1);

%% Components of the binary parameter file "TFG.ForPar"

% AbsTol: Absolute error tolerance for numerical integration.
%         Data type: double, length: 1, units: none
% RelTol: Relative error tolerance for numerical integration.
%         Data type: double, length: 1, units: none

% NumDensityModel: Total number of density models.
%                  Data type: integer, length: 1, units: none
% Fai1: Starting latitudinal boundaries.
%       Data type: double, length: NumDensityModel, units: degree
% Fai2: Ending latitudinal boundaries.
%       Data type: double, length: NumDensityModel, units: degree
% Lamda1: Starting longitudinal boundaries.
%         Data type: double, length: NumDensityModel, units: degree
% Lamda2: Ending longitudinal boundaries.
%         Data type: double, length: NumDensityModel, units: degree
% R1: Starting radial boundaries.
%     Data type: double, length: NumDensityModel, units: user-defined
% R2: Ending radial boundaries.
%     Data type: double, length: NumDensityModel, units: user-defined
% rho: Density of spherical prisms.
%      Data type: double, length: NumDensityModel, units: user-defined

% NumCoordinates: Total number of computational points.
%                 Data type: integer, length: 1, units: none
% Latitude: Latitudes of the computational points.
%           Data type: double, length: NumCoordinates, units: degree
% Longitude: Longitudes of the computational points.
%            Data type: double, length: NumCoordinates, units: degree
% Radius: Radii of the computational points.
%         Data type: double, length: NumCoordinates, units: user-defined

% Notes on units:
% The units of the program's forward-calculated results depend on the units 
% of the input parameters. For example:

% Case 1:
% - G (Newton's gravitational constant): m^3·kg^(-1)·s^(-2)
% - rho (mass density): kg·m^(-3)
% - R1, R2, and Radius: m
% → Resulting units:
%   - Gravitational potential: m^2·s^(-2)
%   - Gravitational acceleration: m·s^(-2)
%   - Gravitational gradient tensor: s^(-2)
%   - Gravitational curvature: m^(-1)·s^(-2)

% Case 2:
% - G (Newton's gravitational constant): m^3·kg^(-1)·s^(-2)
% - rho (mass density): kg·m^(-3)
% - R1, R2, and Radius: km
% → Resulting units:
%   - Gravitational potential: km^2·s^(-2)
%   - Gravitational acceleration: km·s^(-2)
%   - Gravitational gradient tensor: s^(-2)
%   - Gravitational curvature: km^(-1)·s^(-2)

% Reminder:
% - Multiply your calculations by the gravitational constant G.
% - To convert acceleration to commonly used units:
%   1 m·s^(-2) = 10^5 mGal

%%

FileName='TFG.ForPar';
fid=fopen(FileName,'wb+');

fwrite(fid,AbsTol,'double');
fwrite(fid,RelTol,'double');

fwrite(fid,NumDensityModel,'int32');
fwrite(fid,Fai1,'double');
fwrite(fid,Fai2,'double');
fwrite(fid,Lamda1,'double');
fwrite(fid,Lamda2,'double');
fwrite(fid,R1,'double');
fwrite(fid,R2,'double');
fwrite(fid,rho,'double');

fwrite(fid,NumCoordinates,'int32');
fwrite(fid,Longitude,'double');
fwrite(fid,Latitude,'double');
fwrite(fid,Radius,'double');

fclose(fid);





