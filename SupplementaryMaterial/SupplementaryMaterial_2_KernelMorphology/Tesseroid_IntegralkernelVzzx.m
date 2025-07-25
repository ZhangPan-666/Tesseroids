function [KernelVzzx]=Tesseroid_IntegralkernelVzzx(R2,R1,FaiI,LamdaI,r,fai,lamda)

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

KernelVzzx=3.*cosd(FaiI).*Tx.*(R2.^4.*(1-R2.*cos(Phi)).*l1pow5-R1.^4.*(1-R1.*cos(Phi)).*l2pow5)./(l2pow5.*l1pow5);

Index1=Phi<(pi/2-acos(1e-5));
Index2=Phi>(pi/2+acos(1e-5));
KernelVzzx(Index1)=0;
KernelVzzx(Index2)=0;

end