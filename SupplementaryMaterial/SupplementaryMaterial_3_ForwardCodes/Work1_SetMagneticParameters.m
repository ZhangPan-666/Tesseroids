% This Matlab script shows how to generate the parameter files "TFM.ForPar"
% needed for the three cpp main functions:
% Work1_MagneticEstimate.cpp, Work1_MagneticEstimateOpenMP.cpp, and Work1_MagneticEstimateMPI.cpp.
% For more details, see Zhang P. et al.,
% "An alternative approach for accurately calculating gravitational and magnetic fields of a spherical prism."

% For the components of the binary parameter file, please refer to line 47

%%

clc;
clear;
close all;

%%

Fai1=[30;30;40;40];
Fai2=[40;40;50;50];
Lamda1=[110;120;110;120];
Lamda2=[120;130;120;130];

NumPermeabilityModel=size(Fai1,1);

RE=6371.2;

R1=(RE-30).*ones(NumPermeabilityModel,1);
R2=RE.*ones(NumPermeabilityModel,1);

Mx=1.*ones(NumPermeabilityModel,1);
My=1.*ones(NumPermeabilityModel,1);
Mz=1.*ones(NumPermeabilityModel,1);

%%

AbsTol=1e-16;
RelTol=1e-16;

%%

[Longitude,Latitude]=meshgrid(100.05:0.1:140.05,20.05:0.1:60.05);
Longitude=Longitude(:);
Latitude=Latitude(:);

NumCoordinates=size(Longitude,1);
Radius=(RE+1).*ones(NumCoordinates,1);

%% Components of the binary parameter file "TFM.ForPar"

% AbsTol: Absolute error tolerance for numerical integration.
%         Data type: double, length: 1, units: none
% RelTol: Relative error tolerance for numerical integration.
%         Data type: double, length: 1, units: none

% NumPermeabilityModel: Total number of permeability models.
%                       Data type: integer, length: 1, units: none
% Fai1: Starting latitudinal boundaries.
%       Data type: double, length: NumPermeabilityModel, units: degree
% Fai2: Ending latitudinal boundaries.
%       Data type: double, length: NumPermeabilityModel, units: degree
% Lamda1: Starting longitudinal boundaries.
%         Data type: double, length: NumPermeabilityModel, units: degree
% Lamda2: Ending longitudinal boundaries.
%         Data type: double, length: NumPermeabilityModel, units: degree
% R1: Starting radial boundaries.
%     Data type: double, length: NumPermeabilityModel, units: user-defined
% R2: Ending radial boundaries.
%     Data type: double, length: NumPermeabilityModel, units: user-defined
% Mx: Northern component of magnetization intensity of the spherical prisms.
%     Data type: double, length: NumPermeabilityModel, units: user-defined
% My: Eastern component of magnetization intensity of the spherical prisms.
%     Data type: double, length: NumPermeabilityModel, units: user-defined
% Mz: Radial component of magnetization intensity of the spherical prisms.
%     Data type: double, length: NumPermeabilityModel, units: user-defined

% NumCoordinates: Total number of computational points.
%                 Data type: integer, length: 1, units: none
% Latitude: Latitudes of the computational points.
%           Data type: double, length: NumCoordinates, units: degree
% Longitude: Longitudes of the computational points.
%            Data type: double, length: NumCoordinates, units: degree
% Radius: Radii of the computational points.
%         Data type: double, length: NumCoordinates, units: user-defined

% Notes on units:
% The units of the program's forward-calculated magnetic field results
% depend on the units of the input parameters. For example:

% Case 1:
% - mu0 (magnetic permeability of free space): 4·pi·10^(-7) H/m
%   (1 H = 1 m^2·kg·s^(-2)·A^(-2))
% - M (magnetization vector): A/m
% - R1, R2, Radius: m
% → Resulting units:
%   - Magnetic potential: 10^2·m·T
%   - Magnetic vector: 10^2·T
%   - Magnetic gradient tensor: 10^2·m^(-1)·T
%     (1 T = 1 kg·s^(-2)·A^(-1))

% Case 2:
% - mu0: 4·pi·10^(-7) H/m
% - M: A/m
% - R1, R2, Radius: km
% → Resulting units:
%   - Magnetic potential: 10^2·km·T
%   - Magnetic vector: 10^2·T
%   - Magnetic gradient tensor: 10^2·km^(-1)·T
%     (1 T = 1 kg·s^(-2)·A^(-1))

% Unit conversion:
% 1 T = 10^9 nT

%%

FileName='TFM.ForPar';
fid=fopen(FileName,'wb+');

fwrite(fid,AbsTol,'double');
fwrite(fid,RelTol,'double');

fwrite(fid,NumPermeabilityModel,'int32');
fwrite(fid,Fai1,'double');
fwrite(fid,Fai2,'double');
fwrite(fid,Lamda1,'double');
fwrite(fid,Lamda2,'double');
fwrite(fid,R1,'double');
fwrite(fid,R2,'double');
fwrite(fid,Mx,'double');
fwrite(fid,My,'double');
fwrite(fid,Mz,'double');

fwrite(fid,NumCoordinates,'int32');
fwrite(fid,Longitude,'double');
fwrite(fid,Latitude,'double');
fwrite(fid,Radius,'double');

fclose(fid);





