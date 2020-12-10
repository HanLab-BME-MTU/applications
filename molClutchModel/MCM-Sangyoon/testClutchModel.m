%% Experimental values:

% load('experimental data.mat');

%% Common parameters:

nm = 800; %Number of myosin motors, optimal fit 800
fm1 = -2e-12; % Stall force of 1 motor (N)
vu = -110e-9; % Unloaded myosin motor velocity (m/s)
kc = 1; % Clutch spring constant (N/m)
pt = 0.073; % fraction of force experienced by talin 0.073
konv = 1e8; % on-rate of vinculin to unfolded talin
mr = 300*50;  % Maximum integrin density for each integrin
% intadd = 2.4; % Number of integrins added per sq. micron every time reinforcement happens.
a =1700e-9; % Radius of adhesion (m) 1500e-9
ion = 'cm';
ksub = 10.^([-0.1:0.1:2]).*1e-3; %Range of substrate stiffness
kont1 = 2.11e-4; %3.33e-4; % True on-rate (um2/s), 1st integrin type
kont2 = 0; % True on-rate (um2/s), 2nd integrin type
kof1 = 0.9;
kof2 = 1.5;
E = 9*ksub./(4*pi*a);
dint1 = 300; %Density of integrin molecules, type 1 (integrins/um2).
dint2 = 0;   %Density of integrin molecules, type 2 (integrins/um2).
intaddctrl = 24; % At 1000 s: 4 at 100 s: 24
nc10 = 1200; %Number of molecular clutches for 10 ug/ml fn 1200
nc1 = 750; %Number of molecular clutches for 1 ug/ml fn 800
nc100 = 1650; %Number of molecular clutches for 100 ug/ml fn

% 10 ug/ml

% 10 ug/ml depleted

nc = nc10; %Number of molecular clutches
intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(length(ksub),1);
mv = zeros(length(ksub),1);
mnb1 = zeros(length(ksub),1);
mnb2 = zeros(length(ksub),1);
mdint1 = zeros(length(ksub),1);
mdint2 = zeros(length(ksub),1);

%% 
% parfor i=1:length(ksub);
for ii=1:length(ksub)
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/length(ksub)) '% done...'])
end
% mf: Mean force on substrate (N)
% mv: Mean rearward speed (m/s)
% mnb1, mnb2: Mean number of bound clutches for both integrins
% mdint1, mdint2: Mean densities of both integrin types

vdep10 = mv;
mfdep10 = mf;
mdint1dep10 = mdint1;

%% 10 ug/ml - control

nc = nc10; %Number of molecular clutches
intadd = intaddctrl; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(length(ksub),1);
mv = zeros(length(ksub),1);
mnb1 = zeros(length(ksub),1);
mnb2 = zeros(length(ksub),1);
mdint1 = zeros(length(ksub),1);
mdint2 = zeros(length(ksub),1);

parfor ii=1:length(ksub);
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    100*ii./length(ksub)
end;

% matlabpool close;

vctrl10 = mv;
mfctrl10 = mf;
mdint1ctrl10 = mdint1;

Pdep10 = mfdep10./(pi*a.^2);
Pctrl10 = mfctrl10./(pi*a.^2);

figure, 
subplot(2,1,1), semilogx(E,abs(Pctrl10));
hold on, subplot(2,1,1), semilogx(E,abs(Pdep10),'r');
hold on, subplot(2,1,1), errorbar(Eexp(:,1),ctrl10(:,1), ctrl10(:,2),'.b');
hold on, subplot(2,1,1), errorbar(Eexp(:,1),dep10(:,1), dep10(:,2),'.r');
subplot(2,1,2), semilogx(E,1e9*abs(vctrl10));
hold on, subplot(2,1,2), semilogx(E,1e9*abs(vdep10),'r');
hold on, subplot(2,1,2), errorbar(Eexp(:,1),vexpctrl10(:,1), vexpctrl10(:,2),'.b');
hold on, subplot(2,1,2), errorbar(Eexp(:,1),vexpdep10(:,1), vexpdep10(:,2),'.r');
% hold on, 
% subplot(2,1,2),errorbar(Eexp(2:end),vpuroexp, errorvpuro,'.b');
title('10 ug/ml');

figure, 
subplot(2,1,1), semilogx(E,mdint1ctrl10./mdint1ctrl10(1));
hold on, subplot(2,1,1), semilogx(E,mdint1dep10./mdint1dep10(1),'r');
% hold on, subplot(2,1,1), errorbar(Eexp,dint1b1exp, dint1b1exper,'.b');
% subplot(2,1,2), semilogx(E,mdint2b1); 
title('10 ug/ml integrin densities');

% 1 ug/ml

nc = nc1; %Number of molecular clutches best fit: 285
%  1 ug/ml - Talin depleted:


intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(length(ksub),1);
mv = zeros(length(ksub),1);
mnb1 = zeros(length(ksub),1);
mnb2 = zeros(length(ksub),1);
mdint1 = zeros(length(ksub),1);
mdint2 = zeros(length(ksub),1);

% matlabpool open local 3;

parfor ii=1:length(ksub);
    
    [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    100*ii./length(ksub)
end;

vdep1 = mv;
mfdep1 = mf;
mdint1dep1 = mdint1;


% Control - 1 ug/ml

nc = nc1;
intadd = intaddctrl; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(length(ksub),1);
mv = zeros(length(ksub),1);
mnb1 = zeros(length(ksub),1);
mnb2 = zeros(length(ksub),1);
mdint1 = zeros(length(ksub),1);
mdint2 = zeros(length(ksub),1);

parfor ii=1:length(ksub);
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    100*ii./length(ksub)
end;

% matlabpool close;

vctrl1 = mv;
mfctrl1 = mf;
mdint1ctrl1 = mdint1;

Pdep1 = mfdep1./(pi*a.^2);
Pctrl1 = mfctrl1./(pi*a.^2);
figure, 
subplot(2,1,1), semilogx(E,abs(Pctrl1));
hold on, subplot(2,1,1), semilogx(E,abs(Pdep1),'r');
hold on, subplot(2,1,1), errorbar(Eexp(:,1),ctrl1(:,1), ctrl1(:,2),'.b');
hold on, subplot(2,1,1), errorbar(Eexp(:,1),dep1(:,1), dep1(:,2),'.r');
subplot(2,1,2), semilogx(E,1e9*abs(vctrl1));
hold on, subplot(2,1,2), semilogx(E,1e9*abs(vdep1),'r');
% hold on, 
% subplot(2,1,2),errorbar(Eexp(2:end),vpuroexp, errorvpuro,'.b');
title('1 ug/ml');

figure, 
subplot(2,1,1), semilogx(E,mdint1ctrl1./mdint1ctrl1(1));
hold on, subplot(2,1,1), semilogx(E,mdint1dep1./mdint1dep1(1),'r');
% hold on, subplot(2,1,1), errorbar(Eexp,dint1b1exp, dint1b1exper,'.b');
% subplot(2,1,2), semilogx(E,mdint2b1); 
title('1 ug/ml integrin densities');

% 100 ug/ml

nc = nc100; %Number of molecular clutches
%  100 ug/ml - Talin depleted:


intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(length(ksub),1);
mv = zeros(length(ksub),1);
mnb1 = zeros(length(ksub),1);
mnb2 = zeros(length(ksub),1);
mdint1 = zeros(length(ksub),1);
mdint2 = zeros(length(ksub),1);

% matlabpool open local 3;

parfor ii=1:length(ksub);
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    100*ii./length(ksub)
end;

vdep100 = mv;
mfdep100 = mf;
mdint1dep100 = mdint1;


% Control - 100 ug/ml

nc = nc100;
intadd = intaddctrl; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(length(ksub),1);
mv = zeros(length(ksub),1);
mnb1 = zeros(length(ksub),1);
mnb2 = zeros(length(ksub),1);
mdint1 = zeros(length(ksub),1);
mdint2 = zeros(length(ksub),1);

parfor ii=1:length(ksub);
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    100*ii./length(ksub)
end;

matlabpool close;

vctrl100 = mv;
mfctrl100 = mf;
mdint1ctrl100 = mdint1;

Pdep100 = mfdep100./(pi*a.^2);
Pctrl100 = mfctrl100./(pi*a.^2);
figure, 
subplot(2,1,1), semilogx(E,abs(Pctrl100));
hold on, subplot(2,1,1), semilogx(E,abs(Pdep100),'r');
hold on, subplot(2,1,1), errorbar(Eexp(:,1),ctrl100(:,1), ctrl100(:,2),'.b');
hold on, subplot(2,1,1), errorbar(Eexp(:,1),dep100(:,1), dep100(:,2),'.r');
subplot(2,1,2), semilogx(E,1e9*abs(vctrl100));
hold on, subplot(2,1,2), semilogx(E,1e9*abs(vdep100),'r');
% hold on, 
% subplot(2,1,2),errorbar(Eexp(2:end),vpuroexp, errorvpuro,'.b');
title('100 ug/ml');

figure, 
subplot(2,1,1), semilogx(E,mdint1ctrl100./mdint1ctrl100(1));
hold on, subplot(2,1,1), semilogx(E,mdint1dep100./mdint1dep100(1),'r');
% hold on, subplot(2,1,1), errorbar(Eexp,dint1b1exp, dint1b1exper,'.b');
% subplot(2,1,2), semilogx(E,mdint2b1); 
title('100 ug/ml integrin densities');

