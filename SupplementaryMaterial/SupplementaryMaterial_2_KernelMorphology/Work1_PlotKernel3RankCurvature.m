% This Matlab script illustrates the morphology of the surface integration kernel 
% for the 3- rank curvature. 
% For more details, see Zhang P. et al., 
% "An alternative approach for accurately calculating gravitational and magnetic fields of a spherical prism."

% In this script, we consider a single spherical prism with the following coordinate boundaries: 
% radial boundaries r2 = RE, r1 = RE − 30 km; 
% latitudinal boundaries φ2 = 50°, φ1 = 30°; 
% and longitudinal boundaries λ2 = 130° and λ1 = 110°. 
% The computation point is located at φ = 40°, λ = 120°, r = RE + 50 km, i.e., directly above the prism.

% Try plotting the integral kernel with the "mesh" or "contourf" function 
% and observe the morphology of the integral kernel

%%

clc;
clear;
close all;

%%

Fai0=40;
LengthFai=20;
Lamda0=120;
LengthLamda=20;

Fai1=Fai0-LengthFai/2;
Fai2=Fai0+LengthFai/2;
Lamda1=Lamda0-LengthLamda/2;
Lamda2=Lamda0+LengthLamda/2;

RE=6371.2;

R1=RE-30;
R2=RE;

%%

r=RE+50;
fai=40;
lamda=120;

%%

dGrid=0.1;
[Longitude,Latitude]=meshgrid( ...
    Lamda0-LengthLamda/2:dGrid:Lamda0+LengthLamda/2, ...
    Fai0-LengthFai/2:dGrid:Fai0+LengthFai/2);
[Rows,Cols]=size(Longitude);

%%

KernelVxxx=Tesseroid_IntegralkernelVxxx(R2,R1,Latitude,Longitude,r,fai,lamda);

%%

KernelVxxy=Tesseroid_IntegralkernelVxxy(R2,R1,Latitude,Longitude,r,fai,lamda);

%%

KernelVxxz=Tesseroid_IntegralkernelVxxz(R2,R1,Latitude,Longitude,r,fai,lamda);

%%

KernelVxyz=Tesseroid_IntegralkernelVxyz(R2,R1,Latitude,Longitude,r,fai,lamda);

%%

KernelVyyx=Tesseroid_IntegralkernelVyyx(R2,R1,Latitude,Longitude,r,fai,lamda);

%%

KernelVyyy=Tesseroid_IntegralkernelVyyy(R2,R1,Latitude,Longitude,r,fai,lamda);

%%

KernelVyyz=Tesseroid_IntegralkernelVyyz(R2,R1,Latitude,Longitude,r,fai,lamda);

%%

KernelVzzx=Tesseroid_IntegralkernelVzzx(R2,R1,Latitude,Longitude,r,fai,lamda);

%%

KernelVzzy=Tesseroid_IntegralkernelVzzy(R2,R1,Latitude,Longitude,r,fai,lamda);

%%

KernelVzzz=Tesseroid_IntegralkernelVzzz(R2,R1,Latitude,Longitude,r,fai,lamda);