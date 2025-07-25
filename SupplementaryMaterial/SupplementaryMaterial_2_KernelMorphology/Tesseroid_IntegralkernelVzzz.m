function [KernelVzzz]=Tesseroid_IntegralkernelVzzz(R2,R1,FaiI,LamdaI,r,fai,lamda)

h2=r-R2;
hRatio2=h2/r;
h1=r-R1;
hRatio1=h1/r;

R2=1-hRatio2;
R1=1-hRatio1;
Phi=atan2(sqrt((cosd(FaiI).*sind(LamdaI-lamda)).^2+(cosd(fai).*sind(FaiI)-sind(fai).*cosd(FaiI).*cosd(LamdaI-lamda)).^2),...
    sind(fai).*sind(FaiI)+cosd(fai).*cosd(FaiI).*cosd(LamdaI-lamda));

l2pow5=(2.*(2.*sin(Phi/2).^2).*(1-hRatio2)+hRatio2.^2).^(2.5);
l1pow5=(2.*(2.*sin(Phi/2).^2).*(1-hRatio1)+hRatio1.^2).^(2.5);

KernelVzzz=cosd(FaiI).*(...
    -2.*...
    (R2.^3.*l1pow5-R1.^3.*l2pow5)./(l2pow5.*l1pow5)+...
    4.*cos(Phi).*...
    (R2.^4.*l1pow5-R1.^4.*l2pow5)./(l2pow5.*l1pow5)+...
    (1-3.*cos(Phi).^2).*...
    (R2.^5.*l1pow5-R1.^5.*l2pow5)./(l2pow5.*l1pow5));

Index1=Phi<(pi/2-acos(1e-5));
Index2=Phi>(pi/2+acos(1e-5));

if R2<1
    KernelVzzz(Index1)=-2.*cosd(FaiI(Index1)).*(...
        ((1-3.*R2+3.*R2^2)./((1-R2).^3))- ...
        ((1-3.*R1+3.*R1^2)./((1-R1).^3)));
else
    if R1<1
        KernelVzzz(Index1)=2.*cosd(FaiI(Index1))*(...
            ((1-3*R2+3*R2^2)/((1-R2)^3))+ ...
            ((1-3*R1+3*R1^2)/((1-R1)^3)));
    else
        KernelVzzz(Index1)=2.*cosd(FaiI(Index1))*(...
            ((1-3*R2+3*R2^2)/((1-R2)^3))- ...
            ((1-3*R1+3*R1^2)/((1-R1)^3)));
    end
end

KernelVzzz(Index2)=2.*cosd(FaiI(Index2)).*(...
    ((1+3.*R2+3.*R2^2)/((1+R2).^3))- ...
    ((1+3.*R1+3.*R1^2)/((1+R1).^3)));

end