% This Matlab script illustrates the numerical stability of the surface integration kernel 
% for the Vxxz component. 
% For more details, see Zhang P. et al., 
% "An alternative approach for accurately calculating gravitational and magnetic fields of a spherical prism."

% In this numerical experiment, we consider a single spherical prism with the following coordinate boundaries: 
% radial boundaries r2 = RE, r1 = RE − 30 km; 
% latitudinal boundaries φ2 = 50°, φ1 = 30°; 
% and longitudinal boundaries λ2 = 130° and λ1 = 110°. 
% The computation point is located at φ = 40°, λ = 120°, r = RE + 1 km, i.e., directly above the prism.

% To evaluate the kernel behavior, 
% we set the longitude of the running integral point from 119° to 121°,with an interval of 0.001°, 
% while fixing its latitude at 40.001°.

% Figure 1 shows the integration kernel computed using 
% a double-precision numerical environment (in red, showing potential numerical instability) and 
% a symbolic computation environment (in green, representing accurate results). 
% Figure 2 presents the difference between the two.

%%

clc;
clear;
close all;

%% Parameters section No.1
% If you need to make changes, please remember to modify the Parameters section No.2 as well.

RE=6371.2;

R1=RE-30;
R2=RE;

%%

r=RE+1;
fai=40;
lamda=120;

%%

Longitude=119:0.001:121;
Latitude=40.001.*ones(size(Longitude));

%%

FaiI=Latitude;
LamdaI=Longitude;

h2=r-R2;
hRatio2=h2/r;
h1=r-R1;
hRatio1=h1/r;

R2=1-hRatio2;
R1=1-hRatio1;
Phi=atan2(sqrt((cosd(FaiI).*sind(LamdaI-lamda)).^2+(cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda)).^2),...
    sind(fai).*sind(FaiI)+cosd(fai).*cosd(FaiI).*cosd(LamdaI-lamda));
Tx=cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda);

l2pow5=(2.*(2.*sin(Phi/2).^2).*(1-hRatio2)+hRatio2.^2).^(2.5);
l1pow5=(2.*(2.*sin(Phi/2).^2).*(1-hRatio1)+hRatio1.^2).^(2.5);

KernelVxxz1=cosd(FaiI).*(...
    1.*...
    (R2.^3.*l1pow5-R1.^3.*l2pow5)./(l2pow5.*l1pow5)-...
    2.*cos(Phi).*...
    (R2.^4.*l1pow5-R1.^4.*l2pow5)./(l2pow5.*l1pow5)+...
    (1-3.*Tx.^2).*...
    (R2.^5.*l1pow5-R1.^5.*l2pow5)./(l2pow5.*l1pow5));

%% Parameters section No.2

RE=6371.2;

R1=RE-30;
R2=RE;

%%

r=RE+1;
fai=40;
lamda=120;

%%

GuardDigits=32;

FaiI=vpa(Latitude,GuardDigits);
LamdaI=vpa(Longitude,GuardDigits);

h2=vpa(r-R2,GuardDigits);
hRatio2=vpa(h2/r,GuardDigits);
h1=vpa(r-R1,GuardDigits);
hRatio1=vpa(h1/r,GuardDigits);

R2=vpa(1-hRatio2,GuardDigits);
R1=vpa(1-hRatio1,GuardDigits);
Phi=vpa(atan2(sqrt((cosd(FaiI).*sind(LamdaI-lamda)).^2+(cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda)).^2),...
    sind(fai).*sind(FaiI)+cosd(fai).*cosd(FaiI).*cosd(LamdaI-lamda)),GuardDigits);
Tx=vpa(cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda),GuardDigits);

l2pow5=vpa((2.*(2.*sin(Phi/2).^2).*(1-hRatio2)+hRatio2.^2).^(2.5),GuardDigits);
l1pow5=vpa((2.*(2.*sin(Phi/2).^2).*(1-hRatio1)+hRatio1.^2).^(2.5),GuardDigits);

KernelVxxz2=vpa(cosd(FaiI).*(...
    1.*...
    (R2.^3.*l1pow5-R1.^3.*l2pow5)./(l2pow5.*l1pow5)-...
    2.*cos(Phi).*...
    (R2.^4.*l1pow5-R1.^4.*l2pow5)./(l2pow5.*l1pow5)+...
    (1-3.*Tx.^2).*...
    (R2.^5.*l1pow5-R1.^5.*l2pow5)./(l2pow5.*l1pow5)), ...
    GuardDigits);

KernelVxxz2=double(KernelVxxz2);

%%

figure(1);

PP1=plot(Longitude,KernelVxxz1);
PP1.LineWidth=2.0;
PP1.Color=[1,0,0];
PP1.ColorMode='manual';
hold on;

PP1=plot(Longitude,KernelVxxz2);
PP1.LineWidth=2.0;
PP1.Color=[0,1,0];
PP1.ColorMode='manual';
hold off;

AP=gca;
AP.YLim=[-1.1*max(abs(KernelVxxz2)),1.1*max(abs(KernelVxxz2))];

legend('Double-precision computation','Symbolic computation','Location','southeast')

figure(2);

PP1=plot(Longitude,KernelVxxz1-KernelVxxz2);
PP1.LineWidth=2.0;
PP1.Color=[0,0,1];
PP1.ColorMode='manual';
hold off;

legend('Difference','Location','southeast')

