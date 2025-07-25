function [KernelVyy]=Tesseroid_IntegralkernelVyy(R2,R1,FaiI,LamdaI,r,fai,lamda)

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

KernelVyy=cosd(FaiI).*(...
    sin(Alpha).^2.*csc(Phi).^2.*(...
    (-5.*cos(Phi)+3.*cos(Phi).^3).*((l1pow3-l2pow3)./(l2pow3.*l1pow3))+...
    (-3+15.*cos(Phi).^2-6.*cos(Phi).^2.*cos(2.*Phi)).*((R2.*l1pow3-R1.*l2pow3)./(l2pow3.*l1pow3))+...
    (-9.*cos(Phi).^3+3.*cos(Phi).^2.*cos(3.*Phi)).*((R2.^2.*l1pow3-R1.^2.*l2pow3)./(l2pow3.*l1pow3))+...
    (-4+10.*cos(Phi).^2-4.*cos(Phi).^2.*cos(2.*Phi)).*((R2.^3.*l1pow3-R1.^3.*l2pow3)./(l2pow3.*l1pow3)))+...
    cot(Phi).*csc(Phi).*((l1pow1-l2pow1)./(l2pow1.*l1pow1))+...
    (1-cot(Phi).^2).*((R2.*l1pow1-R1.*l2pow1)./(l2pow1.*l1pow1))+...
    (1-3.*Ty.^2).*log((cos(Phi)-R2+l2pow1)./(cos(Phi)-R1+l1pow1)));

Index1=Phi<(pi/2-acos(1e-5));
Index2=Phi>(pi/2+acos(1e-5));

if R2<1
    KernelVyy(Index1)=cosd(FaiI(Index1))*(...
        ((3-4*R2)/(2*(1-R2)^2))- ...
        ((3-4*R1)/(2*(1-R1)^2))+ ...
        log((1-R2)/(1-R1)));
else
    if R1<1
        KernelVyy(Index1)=-cosd(FaiI(Index1))*(...
            ((3-4*R2)/(2*(1-R2)^2))+ ...
            ((3-4*R1)/(2*(1-R1)^2))+ ...
            2*log(r)+ ...
            log((R2-1)*(1-R1)));
    else
        KernelVyy(Index1)=-cosd(FaiI(Index1))*(...
            ((3-4*R2)/(2*(1-R2)^2))- ...
            ((3-4*R1)/(2*(1-R1)^2))+ ...
            log((R2-1)/(R1-1)));
    end
end

KernelVyy(Index2)=-cosd(FaiI(Index2))*(...
    ((3+4*R2)/(2*(1+R2)^2))- ...
    ((3+4*R1)/(2*(1+R1)^2))+ ...
    log((1+R2)/(1+R1)));

end