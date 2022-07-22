clear
clc
nm = 800; %Number of myosin motors, optimal fit 800
fm1 = -2e-12; % Stall force of 1 motor (N)

kc = 5e-4; % Clutch spring constant (N/m)
actinRate=1;

pt = 0.073; % fraction of force experienced by talin 0.073
konv = 1e12; % on-rate of vinculin to unfolded talin
mr = 300*50;  % Maximum integrin density for each integrin
% intadd = 2.4; % Number of integrins added per sq. micron every time reinforcement happens.
a =1700e-9; % Radius of adhesion (m) 1500e-9
ksub = 10.^(-0.1:0.1:1.9).*1e-3; %Range of substrate stiffness
kont1 = 2.11e-4; %3.33e-4; % True on-rate (um2/s), 1st integrin type
kont2 = 0; % True on-rate (um2/s), 2nd integrin type
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
numKsub = length(ksub);
mf = zeros(numKsub,1);
mv = zeros(numKsub,1);
mnb1 = zeros(numKsub,1);
mnb2 = zeros(numKsub,1);
mdint1 = zeros(numKsub,1);
mdint2 = zeros(numKsub,1);
ion = 'mg'; %'cm';
nc = nc10; %Number of molecular clutches
%% testing actin-only mechanosensitivity (blebbi) with no integrin reinforcement
vu = 0; % zero myosin contraction produces zero shortening velocity
v_actin = -12e-9; %-2.6um/min e-6/60 = -4.5e-8 m/s vu = -110e-9; % Unloaded myosin motor velocity (m/s)
intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.

dActin = 1e5; % density of actin at the leading edge #/um

kont1 = 2.11e-3; %increased from 2.11e-4 True on-rate (um2/s), 1st integrin type
kont2 = 0; % True on-rate (um2/s), 2nd integrin type
kof1 = 150;%9; % from 90 previously (5/26/2022)
kof2 = 150; % from 90 previously (5/26/2022)
dint1 = 200; %Density of integrin molecules, type 1 (integrins/um2).
dint2 = 200;   %Density of integrin molecules, type 2 (integrins/um2).

timeTotal = 60;%31*6; % sec
d = 1e-6; % distance from the edge in m.
verbose = 0;
numTrials=5;

Arp_Inh=1;
int_actin=8; 



dActinRange=[dActin];
if Arp_Inh
    dActinRange=[0.5*dActin,0.05*dActin];
    %dActinRange=flip([0:dActin/int_actin:dActin]);
end
v_blebbi_actinSlowdown = zeros(numKsub,numTrials);
mf_blebbi_actinSlowdown = zeros(numKsub,numTrials);
mdint1_blebbi_actinSlowdown = zeros(numKsub,numTrials);
P_blebbi_actinSlowdown = zeros(numKsub,numTrials);

tf=figure;
hold on
legend('Location','bestoutside')
ff=figure;
hold on
legend('Location','bestoutside')


for kk=1:length(dActinRange)
    dActin=dActinRange(kk);
    disp(['Actin Range ' int2str(kk) ' of ' int2str(length(dActinRange))])
    disp(['dActin :: ' int2str(dActin)])

for jj=1:numTrials
    disp(['Starting Trial ' num2str(jj) ' of ' int2str(numTrials)])
    for ii=1:numKsub
        [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
            clutchModelActinElasticityMichels(nm,fm1,vu,nc,dint1,dint2,kont1,...
            kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,dActin,timeTotal,d,verbose,actinRate);
    %     [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
    %        clutchModelNascentAdhesion(nm,fm1,vu,nc,dint1,dint2,kont1,...
    %        kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,dActin);
        mf(ii) = mfi;
        mv(ii) = mvi;
        mnb1(ii) = mnb1i;
        mnb2(ii) = mnb2i;
        mdint1(ii) = mdint1i;
        mdint2(ii) = mdint2i;
        disp(['Trial ' int2str(jj) '::' num2str(100*ii/numKsub) '% done...'])
    end

    v_blebbi_actinSlowdown(:,jj) = mv;
    mf_blebbi_actinSlowdown(:,jj) = mf;
    mdint1_blebbi_actinSlowdown(:,jj) = mdint1;

    P_blebbi_actinSlowdown(:,jj) = mf_blebbi_actinSlowdown(:,jj)/(pi*a^2);
    if verbose
        figure, semilogx(ksub, abs(P_blebbi_actinSlowdown),'o-')
        xlabel('K'), ylabel('Mean traction (Pa)')
        title(['Blebbi, ion: ' ion ', no intadd'])
        %savefig('blebbi_actinElas_Traction.fig')
        figure, semilogx(ksub, abs(v_blebbi_actinSlowdown)*1e9,'o-')
        xlabel('K'), ylabel('Mean flow speed (nm/s)')
        title(['Blebbi, flow speed, ion: ' ion ', no intadd'])
        %savefig('blebbi_actinElas_Flow.fig')
    end
end
figure(tf)
k_basicActin = 1e-6; % basic actin elasiticity: currently totally ambiguous.
errorbar(ksub,abs(mean(P_blebbi_actinSlowdown,2)),(std(P_blebbi_actinSlowdown,0,2))/2,'o-','DisplayName',['k Actin:' num2str(dActinRange(kk)*k_basicActin)]);
set(gca,'XScale','log');
xlabel('K'), ylabel('Mean traction (Pa)')
title(['Blebbi, ion: ' ion ', no intadd', ' Trials:',int2str(numTrials),' Time Period:',int2str(timeTotal)])
figure(ff)
errorbar(ksub,abs(mean(v_blebbi_actinSlowdown,2))*1e6*60,(std(v_blebbi_actinSlowdown,0,2))/2*1e6*60,'o-','DisplayName',['k Actin:' num2str(dActinRange(kk)*k_basicActin)]);
set(gca,'XScale','log');
xlabel('K'), ylabel('Mean flow speed (\mu m/min)')
title(['Blebbi, flow speed, ion: ' ion ', no intadd', ' Trials:',int2str(numTrials),' Time Period:',int2str(timeTotal)])
end