function [KernelVzx]=Tesseroid_IntegralkernelVzx(R2,R1,FaiI,LamdaI,r,fai,lamda)

h2=r-R2;
hRatio2=h2/r;
h1=r-R1;
hRatio1=h1/r;

R2=1-hRatio2;
R1=1-hRatio1;
Phi=atan2(sqrt((cosd(FaiI).*sind(LamdaI-lamda)).^2+(cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda)).^2),...
    sind(fai).*sind(FaiI)+cosd(fai).*cosd(FaiI).*cosd(LamdaI-lamda));
Alpha=atan2(sind(LamdaI-lamda).*cosd(FaiI),cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda));
Tx=cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda);

l2pow1=(2.*(2.*sin(Phi/2).^2).*(1-hRatio2)+hRatio2.^2).^(0.5);
l1pow1=(2.*(2.*sin(Phi/2).^2).*(1-hRatio1)+hRatio1.^2).^(0.5);

l2pow3=(2.*(2.*sin(Phi/2).^2).*(1-hRatio2)+hRatio2.^2).^(1.5);
l1pow3=(2.*(2.*sin(Phi/2).^2).*(1-hRatio1)+hRatio1.^2).^(1.5);

KernelVzx=cosd(FaiI).*(...
    0.5.*cos(Alpha).*(...
    csc(Phi).*(1-3.*cos(2.*Phi)).*((l1pow3-l2pow3)./(l2pow3.*l1pow3))+...
    csc(Phi).*6.*cos(3.*Phi).*((R2.*l1pow3-R1.*l2pow3)./(l2pow3.*l1pow3))+...
    csc(Phi).*3.*(1-2.*cos(2.*Phi)-cos(4.*Phi)).*((R2.^2.*l1pow3-R1.^2.*l2pow3)./(l2pow3.*l1pow3))+...
    (2.*(2.*cos(3.*Phi).*csc(Phi)-cot(Phi))).*((R2.^3.*l1pow3-R1.^3.*l2pow3)./(l2pow3.*l1pow3)))-...
    3.*cos(Phi).*Tx.*log((cos(Phi)-R2+l2pow1)./(cos(Phi)-R1+l1pow1)));

Index1=Phi<(pi/2-acos(1e-5));
Index2=Phi>(pi/2+acos(1e-5));
KernelVzx(Index1)=0;
KernelVzx(Index2)=0;

end