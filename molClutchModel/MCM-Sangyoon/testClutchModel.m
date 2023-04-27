%% Common parameters:
clear
nm = 800; %Number of myosin motors, optimal fit 800
fm1 = -2e-12; % Stall force of 1 motor (N)
vu = -110e-9; % Unloaded myosin motor velocity (m/s)
kc = 1; % Clutch spring constant (N/m)
pt = 0.073; % fraction of force experienced by talin 0.073
konv = 1e8; % on-rate of vinculin to unfolded talin
mr = 300*50;  % Maximum integrin density for each integrin
% intadd = 2.4; % Number of integrins added per sq. micron every time reinforcement happens.
a =1700e-9; % Radius of adhesion (m) 1500e-9
ksub = 0.0025:0.0025:0.05; %10.^(-0.1:0.1:1.5).*1e-3; %Range of substrate stiffness
% ksub = 10.^(-0.1:0.1:2).*1e-3; %Range of substrate stiffness
kont2 = 0; % True on-rate (um2/s), 2nd integrin type
kof2 = 1.5;
E = 9*ksub./(4*pi*a);
dint1 = 300; %Density of integrin molecules, type 1 (integrins/um2).
dint2 = 0;   %Density of integrin molecules, type 2 (integrins/um2).
intaddctrl = 24; % At 1000 s: 4 at 100 s: 24
nc10 = 1200; %Number of molecular clutches for 10 ug/ml fn 1200
nc1 = 750; %Number of molecular clutches for 1 ug/ml fn 800
nc100 = 1650; %Number of molecular clutches for 100 ug/ml fn
tTotal = 100;
% 10 ug/ml

% 10 ug/ml depleted
numKsub = length(ksub);
%% control: intaddctrl = 24; % At 1000 s: 4 at 100 s: 24
close all
ion = 'cm';
nm = 1200; %Number of myosin motors, optimal fit 800
vu = -15e-9; % Unloaded myosin motor velocity (m/s)
nc = nc10; %Number of molecular clutches
% nc = nc1; %Number of molecular clutches
intadd = intaddctrl/10; % Number of integrins added per sq. micron every time reinforcement happens.
v_actin = 0; %m/s
dActin = 5e5; % density of actin at the leading edge #/um
tTotal = 15; % 5 sec to not allow the saturation to a stall force
kont1 = 2.83e-4; % 2.11e-4; % True on-rate (um2/s), 1st integrin type
kof1 = 0.2;
kof2 = 2.5;
dint1 = 300; %Density of integrin molecules, type 1 (integrins/um2).
dint2 = 0;   %Density of integrin molecules, type 2 (integrins/um2).

mf = zeros(numKsub,1);
mv = zeros(numKsub,1);
mnb1 = zeros(numKsub,1);
mnb2 = zeros(numKsub,1);
mdint1 = zeros(numKsub,1);
mdint2 = zeros(numKsub,1);
nExp = 5;
mfGroup = cell(1,nExp);
mvGroup = cell(1,nExp);
mnb1Group = cell(1,nExp);
mnb2Group = cell(1,nExp);
mdint1Group = cell(1,nExp);
mdint2Group = cell(1,nExp);

%Will do five replicate experiments
for p=1:nExp
    for ii=1:numKsub

       [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
           clutchModelNascentAdhesion(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,...
           kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,dActin,tTotal);
        mf(ii) = mfi;
        mv(ii) = mvi;
        mnb1(ii) = mnb1i;
        mnb2(ii) = mnb2i;
        mdint1(ii) = mdint1i;
        mdint2(ii) = mdint2i;
        disp([num2str(100*ii/numKsub) '% done...'])
    end
    mfGroup{p} = mf;
    mvGroup{p} = mv;
    mnb1Group{p} = mnb1;
    mnb2Group{p} = mnb2;
    mdint1Group{p} = mdint1;
    mdint2Group{p} = mdint2;
end
%% plotting traction WT
EkPa = E*1e-3;
vctrl10 = mean(cell2mat(mvGroup),2);
mfctrl10 = mean(cell2mat(mfGroup),2);

vctrl10std = std(cell2mat(mvGroup),0,2);
mfctrl10std = std(cell2mat(mfGroup),0,2);

Pctrl10 = mfctrl10./(pi*a.^2);
Pctrl10err = mfctrl10std./(pi*a.^2);
figure, 
errorbar(EkPa,abs(Pctrl10), Pctrl10err,'ko')

xlabel('E (kPa)'), ylabel('Mean traction (Pa)')
title(['Traction, Control, ion: ' ion ', nc: ' num2str(nc)])

% Fit line
regTract=fit((EkPa)',abs(mean(Pctrl10,2)),'power2');
lim=ylim;
vals=coeffvalues(regTract);
hold on
regT=plot(regTract);
set(regT,'DisplayName',sprintf('T(x)=%.3f*E^{%.3f}+%.3f',vals))
ylim([0,lim(2)]);

% savefig('Control_Traction.fig')
%% plotting flow WT
Vctrl_um = abs(vctrl10)*1e6*60;
Vctrl_um_err = vctrl10std*1e6*60;
figure, %plot(E, abs(vctrl10)*1e6*60,'o')
errorbar(EkPa,Vctrl_um, Vctrl_um_err,"o")
hold on
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,max(Vctrl_um)],...
               'StartPoint',[0.5 1 1]);
md=fittype('a*exp(-b*E)+c','indep','E','options',fo);
% ft = fittype('a*(x-b)^n','problem','n','options',fo);

regSpeed=fit((EkPa'),Vctrl_um,md);
limActin=ylim;
vals=coeffvalues(regSpeed);
regS=plot(regSpeed);
set(regS,'DisplayName',sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',vals))
ylim([0,limActin(2)]);
xlim([0 EkPa(end)])
xlabel('E (kPa)'), ylabel('Mean flow speed (\mum/min)')
title(['Flow, Control, ion: ' ion ', nc: ' num2str(nc)])
% savefig('Control_Flow.fig')

%% testing Roca-cusach's own blebbi simulation via nm
% nm = 180; %Number of myosin motors, optimal fit 800
% vu = -110e-9; % Unloaded myosin motor velocity (m/s)
% v_actin = 0; %-2.6e-6/60; %um/min
% intadd = intaddctrl; % Number of integrins added per sq. micron every time reinforcement happens.
% 
% for ii=1:numKsub
%     
%    [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
%        clutchModelNascentAdhesion(nm,fm1,vu,nc,dint1,dint2,kont1,...
%        kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,1);
%     mf(ii) = mfi;
%     mv(ii) = mvi;
%     mnb1(ii) = mnb1i;
%     mnb2(ii) = mnb2i;
%     mdint1(ii) = mdint1i;
%     mdint2(ii) = mdint2i;
%     disp([num2str(100*ii/numKsub) '% done...'])
% end
% 
% v_blebbiRC = mv;
% mf_blebbiRC = mf;
% mdint1_blebbiRC = mdint1;
% 
% P_blebbiRC = mf_blebbiRC/(pi*a^2);
% figure, semilogx(ksub, abs(P_blebbiRC),'o-')
% xlabel('K'), ylabel('Mean traction (Pa)')
% title('Blebbi, ion: cm, intaddctrl')
% figure, semilogx(ksub, abs(v_blebbiRC),'o-')
% xlabel('K'), ylabel('Mean flow speed')
% title('Blebbi, flow speed, ion: cm, intaddctrl')

%% testing actin-only mechanosensitivity (blebbi) with no integrin reinforcement
% close all
nm=5; vu = 0;
v_actin = -6e-9; %-2.6um/min e-6/60 = -4.5e-8 m/s vu = -110e-9; % Unloaded actin poly velocity (m/s)
intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.
dActin = 7e5; % density of actin at the leading edge #/um. it was 9e5 by 12/04/22
kont1 = 2.11e-4; %increased from 2.11e-4 True on-rate (um2/s), 1st integrin type
kont2 = 0; % True on-rate (um2/s), 2nd integrin type
% kof1 = 1.0; % until 2022/12/03
kof1 = 0.8; % from 90 previously (5/26/2022)
kof2 = 1; % from 90 previously (5/26/2022)
dint1 = 60; %Density of integrin molecules, type 1 (integrins/um2).
dint2 = 0;   %Density of integrin molecules, type 2 (integrins/um2).
ion = 'mg'; %'mg'; %'mg'; % 'cm' doesn't makes sense. Why koff goes up with less force?
tTotal = 10; % sec
d = 1e-6; % distance from the edge in m.
verbose = 1;
mfGroupBBS = cell(1,nExp);
mvGroupBBS = cell(1,nExp);
mnb1GroupBBS = cell(1,nExp);
mnb2GroupBBS = cell(1,nExp);
mdint1GroupBBS = cell(1,nExp);
mdint2GroupBBS = cell(1,nExp);

for p=1:nExp
    for ii=1:numKsub
    %     [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
    %         clutchModelActinElasticity(nm,fm1,vu,nc,dint1,dint2,kont1,...
    %         kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,dActin,timeTotal,d,verbose);
        [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
           clutchModelNascentAdhesion(nm,fm1,vu,nc,dint1,dint2,kont1,...
           kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,dActin,tTotal);
        mf(ii) = mfi;
        mv(ii) = mvi;
        mnb1(ii) = mnb1i;
        mnb2(ii) = mnb2i;
        mdint1(ii) = mdint1i;
        mdint2(ii) = mdint2i;
        disp([num2str(100*ii/numKsub) '% done...'])
    end
    mfGroupBBS{p} = mf;
    mvGroupBBS{p} = mv;
    mnb1GroupBBS{p} = mnb1;
    mnb2GroupBBS{p} = mnb2;
    mdint1GroupBBS{p} = mdint1;
    mdint2GroupBBS{p} = mdint2;
end
%% plot traction bbs vs wt
v_blebbi_actinSlowdown = mean(cell2mat(mvGroupBBS),2); %mv;
mf_blebbi_actinSlowdown = mean(cell2mat(mfGroupBBS),2); %mf;
% mdint1_blebbi_actinSlowdown = mdint1;
mf_blebbi_actinSlowdownErr = std(cell2mat(mfGroupBBS),0,2);

P_blebbi_actinSlowdown = mf_blebbi_actinSlowdown/(pi*a^2);
P_blebbi_actinSlowdownErr = mf_blebbi_actinSlowdownErr/(pi*a^2);

figure, 
errorbar(EkPa, abs(Pctrl10),Pctrl10err,'ko'); hold on


errorbar(EkPa, abs(P_blebbi_actinSlowdown),P_blebbi_actinSlowdownErr,...
    'o','Color',[254/255, 110/255,0])

% Fit line for WT
foWT = fitoptions('Method','NonlinearLeastSquares',...
               'StartPoint',[100 0.3]);
ftTrac=fittype('a*E^(b)-10','indep','E','options',foWT);
regTractWT=fit((EkPa)',abs(Pctrl10),ftTrac); %'power2');
valsWT=coeffvalues(regTractWT);
regTWT=plot(regTractWT,'k');
% regT=plot(regTract,'k');
% Fit line for BBS
foBBS = fitoptions('Method','NonlinearLeastSquares',...
               'StartPoint',[50 0.1]);

ftTrac=fittype('a*E^(b)-10','indep','E','options',foBBS);
regTractBBS=fit((EkPa)',abs(P_blebbi_actinSlowdown),ftTrac); %'power2');
valsBBS=coeffvalues(regTractBBS);
regTBBS=plot(regTractBBS);
regTBBS.Color = [254/255, 110/255,0];
ylim([0,max(abs(Pctrl10))+max(abs(Pctrl10err))]);

% semilogx(ksub, abs(P_blebbi_actinSlowdown),'o','Color',[254/255, 110/255,0])
% xlabel('Spring constant (N/m)'), ylabel('Mean traction (Pa)')
xlabel('Elastic modulus (kPa)'), ylabel('Mean traction (Pa)')
title('Traction WT vs BBS')
legend('WT','BBS',...
    sprintf('T(E)=%.1f*e^{%.2fE}-10',valsWT),...
    sprintf('T(x)=%.1f*E^{%.2fE}-10',valsBBS),'location','best')
xlim([0 20])
% title(['Blebbi, ion: ' ion ', no intadd'])
savefig('./Fig4SimWithErrors/Fig4b.fig')
%% plot for flow WT vs BBS
Vctrl_um = abs(vctrl10)*1e6*60;
Vctrl_um_err = vctrl10std*1e6*60;
VBBS_um = abs(v_blebbi_actinSlowdown)*1e6*60;
v_blebbi_actinSlowdownErr = std(cell2mat(mvGroupBBS),0,2);

VBBS_um_err = v_blebbi_actinSlowdownErr*1e6*60;
figure, 
errorbar(EkPa, abs(Vctrl_um),Vctrl_um_err,'ko'); hold on


errorbar(EkPa, abs(VBBS_um),VBBS_um_err,...
    'o','Color',[254/255, 110/255,0])

% semilogx(ksub, abs(vctrl10)*1e6*60,'ko-'); hold on
% semilogx(ksub, abs(v_blebbi_actinSlowdown)*1e6*60,'o-','Color',[254/255, 110/255,0])
% xlabel('Spring constant (N/m)'), ylabel('Mean flow speed (\mum/min)')
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,max(Vctrl_um)],...
               'StartPoint',[0.5 1 1]);
md=fittype('a*exp(-b*E)+c','indep','E','options',fo);
% ft = fittype('a*(x-b)^n','problem','n','options',fo);

regSpeed=fit((EkPa'),Vctrl_um,md);
limActin=ylim;
valSpWT=coeffvalues(regSpeed);
regS=plot(regSpeed,'k');
set(regS,'DisplayName',sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',valSpWT))
ylim([0,limActin(2)]);
xlim([0 EkPa(end)])

regSpeedBBS=fit((EkPa'),VBBS_um,md);
valsSpBBS=coeffvalues(regSpeedBBS);
regSBBS=plot(regSpeedBBS);
regSBBS.Color = [254/255, 110/255,0];

set(regSBBS,'DisplayName',sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',valsSpBBS))

xlabel('E (kPa)'), ylabel('Mean flow speed (\mum/min)')

title('Flow speed WT vs BBS')
legend('WT','BBS',...
    sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',valSpWT),...
    sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',valsSpBBS),'location','best')
xlim([0 20])
savefig('./Fig4SimWithErrors/Fig4c.fig')
%% title(['Blebbi, flow speed, ion: ' ion ', no intadd'])
save('./Fig4SimWithErrors/Fig4bc.mat')
%% Arp2/3 + myosin inhibition
%Here we keep most parameters and change only F-v sensitivity
v_actin = -6e-9; % same so far
dActin = 9e4; % the only parameter decreased
kof1 = 0.5; % 
tTotal=10;

[mf_ck666,v_ck666,mf_ck666Err,v_ck666Err] = ...
    runClutchModelRangeStiffness(nExp,ksub,nm,fm1,vu,nc,dint1,dint2,kont1,...
           kont2,kof1,kof2,kc,konv,pt,mr,intadd,ion,v_actin,dActin,tTotal);
%% Arp2/3 inhibition: with more unclutches
close all
ion = 'mg';%'mg_earlyslip'; %'cm';
kof1 = 6; % from 30
% tTotal=1000;

[mf_ck666_2,v_ck666_2,mf_ck666_2Err,v_ck666_2Err] = ...
    runClutchModelRangeStiffness(nExp,ksub,nm,fm1,vu,nc,dint1,dint2,kont1,...
           kont2,kof1,kof2,kc,konv,pt,mr,intadd,ion,v_actin,dActin,tTotal);
%% Force
% v_ck666_2 = mv;
% mf_ck666_2 = mf;
% mdint1_ck666_2 = mdint1;

P_ck666 = mf_ck666/(pi*a^2);
P_ck666Err = mf_ck666Err/(pi*a^2);
P_ck666_2 = mf_ck666_2/(pi*a^2);
P_ck666_2Err = mf_ck666_2Err/(pi*a^2);

figure, 
% semilogx(ksub, abs(P_blebbi_actinSlowdown),'o-','Color',[254/255, 110/255,0])
errorbar(EkPa, abs(P_blebbi_actinSlowdown),P_blebbi_actinSlowdownErr,'o','Color',[254/255, 110/255,0])
hold on
regTBBS=plot(regTractBBS);
regTBBS.Color = [254/255, 110/255,0];

% semilogx(ksub, abs(P_ck666),'o-','Color',[0/255, 102/255,204/255])
% P_ck666Err = 0.1*P_ck666;
errorbar(EkPa, abs(P_ck666),P_ck666Err,'o','Color',[0/255, 102/255,204/255])
regTractCK=fit((EkPa)',abs(P_ck666),ftTrac); %'power2');
valsCK=coeffvalues(regTractCK);
regTCK=plot(regTractCK);
regTCK.Color = [0/255, 102/255,204/255];

% semilogx(ksub, abs(P_ck666_2),'o-','Color',[153/255, 0/255,153/255])
% P_ck666_2Err = 0.1*P_ck666_2;
errorbar(EkPa, abs(P_ck666_2),P_ck666_2Err,'o','Color',[153/255, 0/255,153/255])
% xlabel('K'), ylabel('Mean traction (Pa)')
regTractCK2=fit((EkPa)',abs(P_ck666_2),ftTrac); %'power2');
valsCK2=coeffvalues(regTractCK2);
regTCK2=plot(regTractCK2);
regTCK2.Color = [153/255, 0/255,153/255];

xlabel('E (kPa)'), ylabel('Mean traction (Pa)')
% title(['CK666+BBS, ion: ' ion ', no intadd'])
title('CK666+BBS vs BBS')
legend('BBS','CK666+BBS','CK666+BBS+slip','location','best')
legend('BBS', sprintf('T(E)=%.1f*e^{%.2fE}-10',valsBBS),...
    'CK666+BBS', sprintf('T(E)=%.1f*e^{%.2fE}-10',valsCK),...
    'CK666+BBS+slip',sprintf('T(E)=%.1f*e^{%.2fE}-10',valsCK2),'location','best')
savefig('./Fig4SimWithErrors/CK666-BBS-force.fig')
%% flow
VBBS_um = abs(v_blebbi_actinSlowdown)*1e6*60;
VBBS_um_err = v_blebbi_actinSlowdownErr*1e6*60;
figure, 
% semilogx(ksub, abs(v_blebbi_actinSlowdown)*1e6*60,'o-','Color',[254/255, 110/255,0])
errorbar(EkPa, VBBS_um,VBBS_um_err,'o','Color',[254/255, 110/255,0])

hold on
md=fittype('a*exp(-b*E)+c','indep','E','options',fo);
% ft = fittype('a*(x-b)^n','problem','n','options',fo);

regSpeed=fit((EkPa'),VBBS_um,md);
limActin=ylim;
valsBBS=coeffvalues(regSpeed);
regS=plot(regSpeed);
regS.Color = [254/255, 110/255,0];
set(regS,'DisplayName',sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',valsBBS))
ylim([0,limActin(2)]);
xlim([0 EkPa(end)])

v_ck666_um = abs(v_ck666)*1e6*60; 
v_ck666Err_um = abs(v_ck666Err)*1e6*60; 
errorbar(EkPa, v_ck666_um, v_ck666Err_um,'o','Color',[0/255, 102/255,204/255])
regSpeed=fit((EkPa'),v_ck666_um,md);
valsCK1=coeffvalues(regSpeed);
regCK=plot(regSpeed);
regCK.Color = [0/255, 102/255,204/255];
set(regCK,'DisplayName',sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',valsCK1))

% semilogx(ksub, abs(v_ck666)*1e6*60,'o-','Color',[0/255, 102/255,204/255])
v_ck666_2_um = abs(v_ck666_2)*1e6*60; 
v_ck666_2Err_um = abs(v_ck666_2Err)*1e6*60; 
errorbar(EkPa, v_ck666_2_um, v_ck666_2Err_um,'o','Color',[153/255, 0/255,153/255])
regSpeed=fit((EkPa'),v_ck666_2_um,md);
valsCK2=coeffvalues(regSpeed);
regCK2=plot(regSpeed);
regCK2.Color = [153/255, 0/255,153/255];
set(regCK2,'DisplayName',sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',valsCK2))

% semilogx(ksub, abs(v_ck666_2)*1e6*60,'o-','Color',[153/255, 0/255,153/255])
xlabel('E (kPa)'), ylabel('Mean flow speed (\mum/min)')
% xlabel('Spring constant (N/m)'), ylabel('Mean flow speed (\mum/min)')
title('Flow speed BBS vs CK666')
%title(['CK666+BBS, flow speed, ion: ' ion ', no intadd'])
legend('BBS',sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',vals),...
    'CK666+BBS',sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',valsCK1),...
    'CK666+BBS+slip',sprintf('V(x)=%.2f*e^{-%.2f*E}+%.2f',valsBBS),...
    'location','best')
savefig('./Fig4SimWithErrors/CK666-BBS-flow.fig')
save("./Fig4SimWithErrors/Fig4ef.mat")
%% Arp2/3 inhibition: only kon koff control
v_actin = -1e-9; % same so far
kont1 = 2.11e-4; % same as 2.11e-4 True on-rate (um2/s), 1st integrin type
kof1 = 0.5e0; % from 30
% tTotal=1000;
for ii=1:numKsub
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
       clutchModelNascentAdhesion(nm,fm1,vu,nc,dint1,dint2,kont1,...
       kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,dActin);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/numKsub) '% done...'])
end

v_ck666_2 = mv;
mf_ck666_2 = mf;
mdint1_ck666_2 = mdint1;

P_ck666_2 = mf_ck666_2/(pi*a^2);
%% testing actin-elasticity-based model
close all
nc = nc10; %Number of molecular clutches
nm=5; vu = 0;
v_actin = -6e-9; %-2.6um/min e-6/60 = -4.5e-8 m/s vu = -110e-9; % Unloaded actin poly velocity (m/s)
intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.
dActin = 1e11; % density of actin at the leading edge #/um
kont1 = 2.11e-4; %increased from 2.11e-4 True on-rate (um2/s), 1st integrin type
kont2 = 0; % True on-rate (um2/s), 2nd integrin type
kof1 = 3; % from 90 previously (5/26/2022)
kof2 = 1; % from 90 previously (5/26/2022)
dint1 = 100; %Density of integrin molecules, type 1 (integrins/um2).
dint2 = 0;   %Density of integrin molecules, type 2 (integrins/um2).
ion = 'mg'; %'mg'; %'mg'; % 'cm' doesn't makes sense. Why koff goes up with less force?
timeTotal = 10; % milisec
L = 3.2e-11; % m, the length of each actin monomer spring segment. 
d = 1e-6 ; % distance from the edge in m.
eta = 0.1; %viscosity of the dashpot of the F-actin
verbose = 1;
tTotal=10;
folderName=replace(char(datetime(now,'ConvertFrom','datenum')),[" ";":"],'-');
mkdir(folderName);
cd(folderName);

for ii=1:numKsub
    [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
        clutchModelActinElasticity(nm,fm1,nc,dint1,dint2,kont1,...
        kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,dActin,L,...
        timeTotal,d,eta, verbose);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/numKsub) '% done...'])
end

v_blebbi_actinSlowdown = mv;
mf_blebbi_actinSlowdown = mf;
mdint1_blebbi_actinSlowdown = mdint1;

P_blebbi_actinSlowdown = mf_blebbi_actinSlowdown/(pi*a^2);

f1=figure; 
f1.Units = 'inches';
f1.Position(3:4) = [4 3];
% semilogx(ksub, abs(Pctrl10),'ko-'); hold on
% semilogx(ksub, abs(P_blebbi_actinSlowdown),'o-','Color',[254/255, 110/255,0])
plot(ksub, abs(Pctrl10),'ko-'); hold on
plot(ksub, abs(P_blebbi_actinSlowdown),'o-','Color',[254/255, 110/255,0])
xlabel('Spring constant (N/m)'), ylabel('Mean traction (Pa)')
title('Traction WT vs BBS')
legend('WT','BBS','location','best')
% title(['Blebbi, ion: ' ion ', no intadd'])
savefig('blebbi_actinElas_Traction.fig')

f2 = figure;
f2.Units = 'inches';
f2.Position(3:4) = [4 3];
plot(ksub, abs(vctrl10)*1e6*60,'ko-'); hold on
plot(ksub, abs(v_blebbi_actinSlowdown)*1e6*60,'o-','Color',[254/255, 110/255,0])
% semilogx(ksub, abs(vctrl10)*1e6*60,'ko-'); hold on
% semilogx(ksub, abs(v_blebbi_actinSlowdown)*1e6*60,'o-','Color',[254/255, 110/255,0])
xlabel('Spring constant (N/m)'), ylabel('Mean flow speed (\mum/min)')
title('Flow speed WT vs BBS')
legend('WT','BBS','location','best')
savefig('blebbi_actinElas_Flow.fig')
save('workspace.mat')
%% Formin inhibition
vu = 0; % zero myosin contraction produces zero shortening velocity
v_actin = -1e-9; % same so far
dActin = 1e5; % the only parameter decreased
dint1 = 100; %Density of integrin molecules, type 1 (integrins/um2).
ion = 'mg_earlyslip'; % 'cm' doesn't makes sense. Why koff goes up with less force?
kont1 = 5.51e-4; %increased from 2.11e-4 True on-rate (um2/s), 1st integrin type
kof1 = 2.1e0; % from 1.3
tTotal=1000;
parfor ii=1:numKsub
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
       clutchModelNascentAdhesion(nm,fm1,vu,nc,dint1,dint2,kont1,...
       kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,dActin);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/numKsub) '% done...'])
end

v_smif = mv;
mf_smif = mf;
mdint1_smif = mdint1;

P_smif = mf_smif/(pi*a^2);
%% plotting
figure(1)
hold off
plot(E,abs(Pctrl10),'r.-') 
hold on
plot(E,abs(Pdep10),'b.-')
plot(E,abs(P_blebbiRC),'g.-') 
plot(E,abs(P_blebbi_actinSlowdown),'ro-','LineWidth',2) 
plot(E,abs(P_ck666),'co-','LineWidth',2) 
plot(E,abs(P_ck666_2),'ko-','LineWidth',2) 

plot(E,abs(P_smif),'bo-','LineWidth',2) 
% semilogx(E,abs(P_ck666),'co-','LineWidth',2) 
legend('Myosin','Myoin, no FA growth','Blebbi, n_m=180',...
    'No myosin, Actin F-v relationship','No myosin, No Arp2/3 (CK666)',...
    'No myosin, No Arp2/3,vactin=2 (CK666)',...
    'No myosin, No formin (SMIFH2)','Location','best')
ylabel('Traction force')
xlabel('Stiffness')
%% plotting velocity
%This one should be reentered.
h=figure(2);
plot(E,abs(vctrl10),'r.-') 
hold on
% plot(E,abs(v_blebbi),'b.-')
plot(E,abs(v_blebbiRC),'g.-') 
plot(E,abs(v_ck666_2),'ko-','LineWidth',2) 
plot(E,abs(v_smif),'bo-','LineWidth',2) 
ylabel('Actin flow speed')
xlabel('Stiffness')
legend('Myosin','Blebbi, n_m=180',...
    'No myosin, No Arp2/3,vactin=2 (CK666)',...
    'No myosin, No formin (SMIFH2)','Location','best')
try
    figPath = '/Volumes/GoogleDrive/My Drive/Documents/Faculty/Projects/StiffnessSensing/Simulation';
    savefig(h,[figPath,'/flowSpeed_vs_stiffness'])
catch
    savefig(h,[pwd,'/flowSpeed_vs_stiffness'])
    disp('Saving to current directory...')
end

%% testing actin-only mechanosensitivity

vu = 0; % zero myosin contraction produces zero shortening velocity
v_actin = -0.6e-6/60; %-2.6e-6/60; %um/min
nc = nc10; %Number of molecular clutches
intadd = intaddctrl; % Number of integrins added per sq. micron every time reinforcement happens.

for ii=1:numKsub
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
       clutchModelNascentAdhesion(nm,fm1,vu,nc,dint1,dint2,kont1,...
       kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/numKsub) '% done...'])
end

v_blebbi = mv;
mf_blebbi = mf;
mdint1_blebbi = mdint1;

P_blebbi = mf_blebbi/(pi*a^2);
%% talin2 shRNA: intadd=0: same as chan and odde
ion = 'cm'; %'cm';
nc = nc10; %Number of molecular clutches
v_actin = -1.5e-9; %-2.6um/min e-6/60 = -4.5e-8 m/s vu = -110e-9; % Unloaded myosin motor velocity (m/s)
intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.
dActin = 1e4; % density of actin at the leading edge #/um
% nm = 750; fm1 = -2e-12; vu = -120e-9; 
% kc = 5e-3; % Clutch spring constant (N/m)
% pt = 0.073; % fraction of force experienced by talin 0.073
% konv = 0; % on-rate of vinculin to unfolded talin
% mr = 300*50;  % Maximum integrin density for each integrin
% a = 1700e-9; % Radius of adhesion (m) 1500e-9
% ion = 'cm';
% kont1 = 0.001; %3.33e-4; % True on-rate (um2/s), 1st integrin type
% kont2 = 0; % True on-rate (um2/s), 2nd integrin type
% kof1 = 0.9;
% kof2 = 0.8;
% ksub = 10.^(-3:0.5:2); %Range of substrate stiffness
% E = 9*ksub./(4*pi*a);
% dint1 = 30; %Density of integrin molecules, type 1 (integrins/um2).
% dint2 = 0;   %Density of integrin molecules, type 2 (integrins/um2).

for ii=1:numKsub
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchModelNascentAdhesion(nm,...
       fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,...
       mr,intadd,ion,v_actin,dActin, tTotal);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/numKsub) '% done...'])
end
% mf: Mean force on substrate (N)
% mv: Mean rearward speed (m/s)
% mnb1, mnb2: Mean number of bound clutches for both integrins
% mdint1, mdint2: Mean densities of both integrin types

vdep10 = mv;
mfdep10 = mf;
mdint1dep10 = mdint1;
Pdep10 = mfdep10./(pi*a.^2);
figure, semilogx(ksub, abs(Pdep10),'o-')
xlabel('K'), ylabel('Mean traction (Pa)')
title('talin2 shRNA, ion: cm')
%% 1 ug/ml

nc = nc1; %Number of molecular clutches best fit: 285
%  1 ug/ml - Talin depleted:


intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(numKsub,1);
mv = zeros(numKsub,1);
mnb1 = zeros(numKsub,1);
mnb2 = zeros(numKsub,1);
mdint1 = zeros(numKsub,1);
mdint2 = zeros(numKsub,1);

% matlabpool open local 3;
parfor ii=1:numKsub
    
    [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/numKsub) '% done...'])
end

vdep1 = mv;
mfdep1 = mf;
mdint1dep1 = mdint1;


%% Control - 1 ug/ml: FN density

nc = nc1;
intadd = intaddctrl; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(numKsub,1);
mv = zeros(numKsub,1);
mnb1 = zeros(numKsub,1);
mnb2 = zeros(numKsub,1);
mdint1 = zeros(numKsub,1);
mdint2 = zeros(numKsub,1);

parfor ii=1:numKsub
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/numKsub) '% done...'])
end

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

%% 100 ug/ml of FN density

nc = nc100; %Number of molecular clutches
%  100 ug/ml - Talin depleted:


intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(numKsub,1);
mv = zeros(numKsub,1);
mnb1 = zeros(numKsub,1);
mnb2 = zeros(numKsub,1);
mdint1 = zeros(numKsub,1);
mdint2 = zeros(numKsub,1);

% matlabpool open local 3;

parfor ii=1:numKsub
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/numKsub) '% done...'])
end

vdep100 = mv;
mfdep100 = mf;
mdint1dep100 = mdint1;

%% Control - 100 ug/ml FN density

nc = nc100;
intadd = intaddctrl; % Number of integrins added per sq. micron every time reinforcement happens.

mf = zeros(numKsub,1);
mv = zeros(numKsub,1);
mnb1 = zeros(numKsub,1);
mnb2 = zeros(numKsub,1);
mdint1 = zeros(numKsub,1);
mdint2 = zeros(numKsub,1);

parfor ii=1:numKsub
    
   [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion);
    mf(ii) = mfi;
    mv(ii) = mvi;
    mnb1(ii) = mnb1i;
    mnb2(ii) = mnb2i;
    mdint1(ii) = mdint1i;
    mdint2(ii) = mdint2i;
    disp([num2str(100*ii/numKsub) '% done...'])
end

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

