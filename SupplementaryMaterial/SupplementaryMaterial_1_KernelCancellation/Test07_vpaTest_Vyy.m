% This Matlab script illustrates the numerical stability of the surface integration kernel 
% for the Vyy component. 
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
Alpha=atan2(sind(LamdaI-lamda).*cosd(FaiI),cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda));
Ty=cosd(FaiI).*sind(LamdaI-lamda);

l2pow1=(2.*(2.*sin(Phi/2).^2).*(1-hRatio2)+hRatio2.^2).^(0.5);
l1pow1=(2.*(2.*sin(Phi/2).^2).*(1-hRatio1)+hRatio1.^2).^(0.5);

l2pow3=(2.*(2.*sin(Phi/2).^2).*(1-hRatio2)+hRatio2.^2).^(1.5);
l1pow3=(2.*(2.*sin(Phi/2).^2).*(1-hRatio1)+hRatio1.^2).^(1.5);

KernelVyy1=cosd(FaiI).*(...
    sin(Alpha).^2.*csc(Phi).^2.*(...
    (-5.*cos(Phi)+3.*cos(Phi).^3).*((l1pow3-l2pow3)./(l2pow3.*l1pow3))+...
    (-3+15.*cos(Phi).^2-6.*cos(Phi).^2.*cos(2.*Phi)).*((R2.*l1pow3-R1.*l2pow3)./(l2pow3.*l1pow3))+...
    (-9.*cos(Phi).^3+3.*cos(Phi).^2.*cos(3.*Phi)).*((R2.^2.*l1pow3-R1.^2.*l2pow3)./(l2pow3.*l1pow3))+...
    (-4+10.*cos(Phi).^2-4.*cos(Phi).^2.*cos(2.*Phi)).*((R2.^3.*l1pow3-R1.^3.*l2pow3)./(l2pow3.*l1pow3)))+...
    cot(Phi).*csc(Phi).*((l1pow1-l2pow1)./(l2pow1.*l1pow1))+...
    (1-cot(Phi).^2).*((R2.*l1pow1-R1.*l2pow1)./(l2pow1.*l1pow1))+...
    (1-3.*Ty.^2).*log((cos(Phi)-R2+l2pow1)./(cos(Phi)-R1+l1pow1)));

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
Alpha=vpa(atan2(sind(LamdaI-lamda).*cosd(FaiI),cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda)),GuardDigits);
Ty=vpa(cosd(FaiI).*sind(LamdaI-lamda),GuardDigits);

l2pow1=vpa((2.*(2.*sin(Phi/2).^2).*(1-hRatio2)+hRatio2.^2).^(0.5),GuardDigits);
l1pow1=vpa((2.*(2.*sin(Phi/2).^2).*(1-hRatio1)+hRatio1.^2).^(0.5),GuardDigits);

l2pow3=vpa((2.*(2.*sin(Phi/2).^2).*(1-hRatio2)+hRatio2.^2).^(1.5),GuardDigits);
l1pow3=vpa((2.*(2.*sin(Phi/2).^2).*(1-hRatio1)+hRatio1.^2).^(1.5),GuardDigits);

KernelVyy2=vpa(cosd(FaiI).*(...
    sin(Alpha).^2.*csc(Phi).^2.*(...
    (-5.*cos(Phi)+3.*cos(Phi).^3).*((l1pow3-l2pow3)./(l2pow3.*l1pow3))+...
    (-3+15.*cos(Phi).^2-6.*cos(Phi).^2.*cos(2.*Phi)).*((R2.*l1pow3-R1.*l2pow3)./(l2pow3.*l1pow3))+...
    (-9.*cos(Phi).^3+3.*cos(Phi).^2.*cos(3.*Phi)).*((R2.^2.*l1pow3-R1.^2.*l2pow3)./(l2pow3.*l1pow3))+...
    (-4+10.*cos(Phi).^2-4.*cos(Phi).^2.*cos(2.*Phi)).*((R2.^3.*l1pow3-R1.^3.*l2pow3)./(l2pow3.*l1pow3)))+...
    cot(Phi).*csc(Phi).*((l1pow1-l2pow1)./(l2pow1.*l1pow1))+...
    (1-cot(Phi).^2).*((R2.*l1pow1-R1.*l2pow1)./(l2pow1.*l1pow1))+...
    (1-3.*Ty.^2).*log((cos(Phi)-R2+l2pow1)./(cos(Phi)-R1+l1pow1))),...
    GuardDigits);

KernelVyy2=double(KernelVyy2);

%%

figure(1);

PP1=plot(Longitude,KernelVyy1);
PP1.LineWidth=2.0;
PP1.Color=[1,0,0];
PP1.ColorMode='manual';
hold on;

PP1=plot(Longitude,KernelVyy2);
PP1.LineWidth=2.0;
PP1.Color=[0,1,0];
PP1.ColorMode='manual';
hold off;

AP=gca;
AP.YLim=[-1.1*max(abs(KernelVyy2)),1.1*max(abs(KernelVyy2))];

legend('Double-precision computation','Symbolic computation','Location','northeast')

figure(2);

PP1=plot(Longitude,KernelVyy1-KernelVyy2);
PP1.LineWidth=2.0;
PP1.Color=[0,0,1];
PP1.ColorMode='manual';
hold off;

legend('Difference','Location','southeast')