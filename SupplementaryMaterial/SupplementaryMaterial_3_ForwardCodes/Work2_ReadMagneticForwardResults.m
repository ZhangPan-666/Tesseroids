% This Matlab script shows how to read the result files "Result_TFM.dat"
% generated from three cpp main functions:
% Work1_MagneticEstimate.cpp, Work1_MagneticEstimateOpenMP.cpp, and Work1_MagneticEstimateMPI.cpp.
% For more details, see Zhang P. et al.,
% "An alternative approach for accurately calculating gravitational and magnetic fields of a spherical prism."

%%

clc;
clear;
close all;

%%

RE=6371.2;

fid=fopen('TFM.ForPar','rb+');

AbsTol=fread(fid,1,'double');
RelTol=fread(fid,1,'double');

NumDensityModel=fread(fid,1,'int32');
Fai1=fread(fid,NumDensityModel,'double');
Fai2=fread(fid,NumDensityModel,'double');
Lamda1=fread(fid,NumDensityModel,'double');
Lamda2=fread(fid,NumDensityModel,'double');
R1=fread(fid,NumDensityModel,'double');
R2=fread(fid,NumDensityModel,'double');
Mx=fread(fid,NumDensityModel,'double');
My=fread(fid,NumDensityModel,'double');
Mz=fread(fid,NumDensityModel,'double');

fclose(fid);

%%

fid=fopen('Result_TFM.dat','rb+');

NumCoordinates=fread(fid,1,'int32');

Longitude=fread(fid,NumCoordinates,'double');
Latitude=fread(fid,NumCoordinates,'double');
Radius=fread(fid,NumCoordinates,'double');

V=fread(fid,NumCoordinates,'double');

Vx=fread(fid,NumCoordinates,'double');
Vy=fread(fid,NumCoordinates,'double');
Vz=fread(fid,NumCoordinates,'double');

Vxx=fread(fid,NumCoordinates,'double');
Vxy=fread(fid,NumCoordinates,'double');
Vyy=fread(fid,NumCoordinates,'double');
Vzx=fread(fid,NumCoordinates,'double');
Vzy=fread(fid,NumCoordinates,'double');
Vzz=fread(fid,NumCoordinates,'double');

fclose(fid);





