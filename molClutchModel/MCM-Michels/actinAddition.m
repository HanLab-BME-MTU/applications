figure()
hold on
dActin = 0.5e5;

d = 2e-6;
L = 2e-9; % m, the length of each actin monomer spring segment. 
Norg = d/L; % The number of actin springs in-between the membrane and adhesion
Nnew = 0; % Newly-added actin springs at the membrane in front of the adhesion
Nall = Norg + Nnew;
kB = 1.38064852e-23; %m2 kg s-2 K-1
T = 278; %K
c = 0.8;
pot=1;
Nn=[];
mA=[];
fA=[];
j=1;
for i=[0:dActin/10:dActin]
    k_actin=i;
    k0=k_actin;
    C_actin = kB*T*c*dActin; %constant for force-velocity relationship in actin: This is assumption for now 
    R = 1e-6; % m, the radius of curvature of edge. Given normal cell, it can be ~ 10-30 um
    Fs_actin = C_actin/(4*R);
    Nnew_cur=0;
    maxActin_0=20;
    maxActin=maxActin_0;
    boundTime=1;

    
    while (Nnew_cur<100) %actin addition
        %FcNext=max(abs(Fc))+Nnew_cur/(Nall+Nnew_cur)*L*k_actin;
        Fa0_l=(Nnew+Nnew_cur)*L*k_actin/(Nall+Nnew_cur);
        Fa_l=Fa0_l*(1-exp(-k0*boundTime/pot));
        maxActin=maxActin_0*(1-Fa_l/Fs_actin);
        Nnew_cur=Nnew_cur+1;

    end
    Nn(j)=Nnew_cur;
    mA(j)=maxActin;
    fA(j)=Fa_l;
    j=j+1;
end


