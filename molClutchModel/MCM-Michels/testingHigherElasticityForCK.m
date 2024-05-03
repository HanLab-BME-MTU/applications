% testingHigherElasticityForCK.m is a script that test higher kactin from
% BBS+CK666 inhibition case using clutchModelActinElasticityMichels
% function. 
% The main parameters that will be modulated will be:
% dActinRange=[11*0.01e11 0.015e11];
% etaRange=[0.175 0.175];
% offset=["0.10" "0.10"];
% kof1Range=[1.7 3];
% NnMaxRange=[1,1];
% The main one is dActinRange. I will change it to 
% dActinRange=[11e9 4e9 4e9];
% With this, I might have to adjust eta and kof
% I'll try with the current set first
% etaRange=[0.175    1.0800   0.8];
% kof1Range=[1.7 3];

%% initialization
clear
close all
clc

figures=[];
figNames=[];

nm = 5; %Number of myosin motors, optimal fit 800
fm1 = -2e-8; % Stall force of 1 motor (N)

kc = 1; % Clutch spring constant (N/m)
actinRate=2;

pt = 0.073; % fraction of force experienced by talin 0.073
konv = 1e8; % on-rate of vinculin to unfolded talin
mr = 300*50;  % Maximum integrin density for each integrin
% intadd = 2.4; % Number of integrins added per sq. micron every time reinforcement happens.
a =1700e-9; % Radius of adhesion (m) 1500e-9

% ksub = 10.^(-0.1:0.1:1.5).*1e-3; %Range of substrate stiffness (full)

ksub=[0.6 1.3 2.6 6 12.7]*(4*pi*a)/9*10^3; %Range of substrate stiffness (match experiments)


intaddctrl = 24; % At 1000 s: 4 at 100 s: 24
nc10 = 1200; %Number of molecular clutches for 10 ug/ml fn 1200
nc1 = 750; %Number of molecular clutches for 1 ug/ml fn 800
nc100 = 1650; %Number of molecular clutches for 100 ug/ml fn
% 10 ug/ml

% 10 ug/ml depleted
%ksub=[1e-3,1e-1];

E = 9*ksub./(4*pi*a);
numKsub = length(ksub);
mf = zeros(numKsub,1);
mv = zeros(numKsub,1);
mnb1 = zeros(numKsub,1);
mnb2 = zeros(numKsub,1);
mdint1 = zeros(numKsub,1);
mdint2 = zeros(numKsub,1);
mfc = zeros(numKsub,1);
ion = 'mg'; %'mg_earlyslip'; %'cm';
nc = nc10; %Number of molecular clutches
%% testing actin-only mechanosensitivity (blebbi) with no integrin reinforcement
vu = 0; % zero myosin contraction produces zero shortening velocity
v_actin = -6e-9; %-2.6um/min e-6/60 = -4.5e-8 m/s vu = -110e-9; % Unloaded myosin motor velocity (m/s)
intadd = 0; % Number of integrins added per sq. micron every time reinforcement happens.

 

kont1 = 2.11e-4; %increased from 2.11e-4 True on-rate (um2/s), 1st integrin type
kont2 = 0; % True on-rate (um2/s), 2nd integrin type
kof1 = 1.5;%0.4; %250;%250;%200;%9; % from 90 previously (5/26/2022)
kof2 = 1;%250;%200; % from 90 previously (5/26/2022)
dint1 = 100; %Density of integrin molecules, type 1 (integrins/um2).
dint2 = 0;   %Density of integrin molecules, type 2 (integrins/um2).
d = 1e-6; % distance from the edge in m.
k_basicActin = 1e-6; % basic actin elasiticity: currently totally ambiguous.



%% TRIAL SETTINGS 

slip=0; %if zero uses previous time steps velocity when clutch slips, prevents jumping velocity
L = 32e-9; %m, thus 32 nm %3.2e-11; %.3e-12; % m, the length of each actin monomer spring segment. 
dActin = 1e11; % density of actin at the leading edge #/um
eta=0.1;%0.0003;   

timeTotal = 200; % sec
verbose = 0; % time based figures for each trial, if enabled will open and save 
showFreq=0; % if 1, will open frequency plots
numTrials=2; % 5 % number of trials to average for figures
NnMax=1;


%% Iterate dActin
dActinIt=0; % if 1 will iterate over the specified parameter
int_actin=15; % number of steps to iterate for 
maxA=1.5; % multiplier for value 

%% Iterate k Clutch
kcIt=0;
int_kc=10;
maxC=2;

%% Iterate eta clutch
etaIt=1;
int_eta=5;
maxP=20;
kcRange=[kc];

%% set arrays for iteration

v_blebbi_actinSlowdown = zeros(numKsub,numTrials);
mf_blebbi_actinSlowdown = zeros(numKsub,numTrials);
mdint1_blebbi_actinSlowdown = zeros(numKsub,numTrials);
P_blebbi_actinSlowdown = zeros(numKsub,numTrials);

tf=figure;
figures=[figures(:);tf];
figNames=[figNames(:);"Mean_Traction"];
hold on
legend('Location','bestoutside')
ff=figure;
figures=[figures(:);ff];
figNames=[figNames(:);"Mean_Flow_Speed"];
hold on
legend('Location','bestoutside')
cc=figure;
figures=[figures(:);cc];
figNames=[figNames(:);"Mean_Clutch_Force"];
hold on
legend('Location','bestoutside')

%% Params for blebbi and ck666

dActinRange=[11e9 4e9 7e9];
etaRange=[0.175 0.7 1.2];

offset=["0.10" "0.10" "0.10"];
kof1Range=[1.7 2 2];
NnMaxRange=[1,1,1];

%% Params for ck666, smifh2, and LatA

% % dActinRange=[0.1500e10    0.115e10    0.1e10];
% dActinRange=[1.5e9    0.5e9    0.3e9];
% etaRange=[0.175    1.0800   0.8];
% 
% offset=["0.1" "0.02" "0.02"];
% kof1Range=[3 3 3];
% NnMaxRange=[1,1,0.5];

%% params for blebbi figs
% 
% dActinRange=[8*0.01e11];
% etaRange=[0.1];
% 
% offset=["0.1" ];
% kof1Range=[1.7 ];
% NnMaxRange=[1];



%% simulation
cm=prism(length(kcRange)*length(dActinRange)*length(etaRange));
cmi=1;
for kk=1:length(dActinRange)
    dActin=dActinRange(kk);
    eta=etaRange(kk);
    kof1=kof1Range(kk);
    NnMax=NnMaxRange(kk);
    disp(['Actin Range ' int2str(kk) ' of ' int2str(length(dActinRange))])
    disp(['dActin :: ' int2str(dActin)])
    
    %eta=0.1-(1e-12*dActin)^2;

    for jj=1:numTrials
        disp(['Starting Trial ' num2str(jj) ' of ' int2str(numTrials)])
        for ii=1:numKsub
            [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i,mfci,sffti{1}] = ...
                clutchModelActinElasticityMichels(nm,fm1,nc,dint1,dint2,kont1,...
                kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,dActin,timeTotal,d,verbose,actinRate,eta,slip,L,NnMax);
        %     [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
        %        clutchModelNascentAdhesion(nm,fm1,vu,nc,dint1,dint2,kont1,...
        %        kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,dActin);
            mf(ii) = mfi;
            mv(ii) = mvi;
            mnb1(ii) = mnb1i;
            mnb2(ii) = mnb2i;
            mdint1(ii) = mdint1i;
            mdint2(ii) = mdint2i;
            mfc(ii)=mfci;
            sfftt(ii)=sffti;
            disp(['Trial ' int2str(jj) '::' num2str(100*ii/numKsub) '% done...'])
        end
    
        v_blebbi_actinSlowdown(:,jj) = mv;
        mf_blebbi_actinSlowdown(:,jj) = mf;
        mdint1_blebbi_actinSlowdown(:,jj) = mdint1;
        sffta(:,jj)=sfftt;

        P_blebbi_actinSlowdown(:,jj) = mf_blebbi_actinSlowdown(:,jj)/(pi*a^2);
        if verbose
    %         figure, semilogx(ksub, abs(P_blebbi_actinSlowdown),'o-')
    %         xlabel('K'), ylabel('Mean traction (Pa)')
    %         title(['Blebbi, ion: ' ion ', no intadd'])
    %         %savefig('blebbi_actinElas_Traction.fig')
    %         figure, semilogx(ksub, abs(v_blebbi_actinSlowdown)*1e9,'o-')
    %         xlabel('K'), ylabel('Mean flow speed (nm/s)')
    %         title(['Blebbi, flow speed, ion: ' ion ', no intadd'])
    %         %savefig('blebbi_actinElas_Flow.fig')
        end
        figure(cc)
        %semilogx(E*10^-3,abs(mfc),'Marker','o','DisplayName',['k Actin:' num2str(dActinRange(kk)*k_basicActin) ' kc:' num2str(kc) ' \eta:' num2str(eta)],'Color',cm(cmi,:))
        plot(E*10^-3,abs(mfc),'Marker','o','DisplayName',['k Actin:' num2str(dActinRange(kk)*k_basicActin) ' kc:' num2str(kc) ' \eta:' num2str(eta)],'Color',cm(cmi,:))
        xlabel('E (kPa)'), ylabel('Force Clutch (N)')
        title(['Mean Force Clutch, ion: ' ion ', no intadd'])
    end


    figure(tf)
    errorbar(E*10^-3,abs(mean(P_blebbi_actinSlowdown,2)),(std(P_blebbi_actinSlowdown,0,2)),'o-','DisplayName',['k Actin:' num2str(dActinRange(kk)*k_basicActin) ' \eta:' num2str(eta)],'Color',cm(cmi,:),'LineStyle','none');
    regTract=fit((E*10^-3)',abs(mean(P_blebbi_actinSlowdown,2)),'power1');
    lim=ylim;
    vals=coeffvalues(regTract);
    regT=plot(regTract);
    set(regT,'DisplayName',sprintf('F(E)=%.3f*E^{%.3f}',vals))
    set(regT,'Color',cm(cmi,:));
    ylim([0,lim(2)]);
    %set(gca,'XScale','log');
    xlabel('E (kPa)'), ylabel('Mean traction (Pa)')
    %title(['Blebbi, ion: ' ion ', no intadd', ' Trials:',int2str(numTrials),' Time Period:',int2str(timeTotal)])
    figure(ff)
    errorbar(E*10^-3,abs(mean(v_blebbi_actinSlowdown,2))*1e6*60,(std(v_blebbi_actinSlowdown,0,2))*1e6*60,'o-','DisplayName',['k Actin:' num2str(dActinRange(kk)*k_basicActin) ' \eta:' num2str(eta)],'Color',cm(cmi,:),'LineStyle','none');
    %set(gca,'XScale','log');
    md=fittype(['a*exp(b*t)+' char(offset(kk)) ],'indep','t');
    regSpeed=fit((E*10^-3)',abs(mean(v_blebbi_actinSlowdown,2))*1e6*60,md);
    lim=ylim;
    vals=coeffvalues(regSpeed);
    regS=plot(regSpeed);
    set(regS,'DisplayName',sprintf(['V(E)=%.3fe^{%.3f*E}+' char(offset(kk))],vals))
    set(regS,'Color',cm(cmi,:));
    ylim([0,lim(2)]);
    xlabel('E (kPa)'), ylabel('Mean flow speed (\mu m/min)')
    %title(['Blebbi, flow speed, ion: ' ion ', no intadd', ' Trials:',int2str(numTrials),' Time Period:',int2str(timeTotal)])

    drawnow
    if (showFreq)
        sfft=cell(length(sffta),2);
        sfft(1:length(sffta))=(sffta{1,1}(1));
        for tt=[1:length(sfft)]
            powers=zeros(numTrials,length(sffta{1,1}{2}));
            for ttt=[1:numTrials]
                powers(ttt,:)=sffta{tt,ttt}{2}(1:end);
            end
            sfft{tt,2}=mean(powers,1);
        end
        numfreq=25;
        spectrum=zeros(numfreq,length(sfft));

        for tt =[1:length(sfft)]
            spectrum(:,tt)=normalize(flip(sfft{tt,2}(1:numfreq)));
        end
        % 2d histogram/density plot
        spec=figure;
        figures=[figures(:);spec];
        figNames=[figNames(:);string(['Freq_Spectrum_' num2str(dActinRange(kk)*k_basicActin) '_' num2str(kc)  '_' num2str(eta)])];
        numX=15;
        numY=numfreq;
        image(spectrum,'CDataMapping','scaled');
        colormap(spec, "turbo")
        cb=colorbar(spec.CurrentAxes,"eastoutside");
        cb.Label.String="Power";
        ylabel("Frequency (Hz)");
        xlabel("E (kPa)");
        flabels=cellfun(@(f) sprintf('%f',f),num2cell(sfft{1,1}(1:numfreq)),'UniformOutput',false);
        elabels=cellfun(@(f) sprintf('%0.1f',f),num2cell(E),'UniformOutput',false);
        elabels=elabels(1:floor(length(E)/numX):end);
        flabels=flabels(1:floor(length(flabels)/numY):end);
        set(spec.CurrentAxes,'Ytick',linspace(0,numfreq,numY));
        set(spec.CurrentAxes,'YTickLabels',flabels(end:-1:1));
        set(spec.CurrentAxes,'XTick',linspace(0,length(E),numX));
        set(spec.CurrentAxes,'XTickLabels',elabels(1:end));
        
        ap=figure;
        figures=[figures(:);ap];
        figNames=[figNames(:);string(['Freq_Area_' num2str(dActinRange(kk)*k_basicActin) '_' num2str(kc)  '_' num2str(eta)])];
        hold on;
        cmap=turbo(length(E));
        for aa=[1:length(E)]
            area(sfft{aa,1}(1:numfreq),normalize(sfft{aa,2}(1:numfreq)),'FaceColor',cmap(aa,:),'FaceAlpha',.5,'DisplayName',[num2str(E(aa))])
        end
        xlabel('Frequency (Hz)');
        ylabel('Normalised Power');
        apl=legend;
        title(apl,['E (kPa)']);
        hold off;
        drawnow;
    end
    cmi=cmi+1;
end
%%save workspace & figs (does not save figures generated within model
folderName=replace(char(datetime(now,'ConvertFrom','datenum')),[" ";":"],'-');
mkdir(folderName);
cd(folderName);
for fig = [1:length(figures)]
    try
        exportgraphics(figures(fig),[char(figNames(fig)) '.jpg'])
        savefig(figures(fig),[char(figNames(fig)) '.fig'])
        print(figures(fig),char(figNames(fig)),'-djpeg')
    catch
        savefig(figures(fig),[char(figNames(fig)) '.fig'])
        print(figures(fig),char(figNames(fig)),'-djpeg')
    end
end
save('workspace.mat')
cd ..

