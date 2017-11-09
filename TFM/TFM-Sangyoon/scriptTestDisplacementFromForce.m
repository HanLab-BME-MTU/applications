% scriptTestDisplacementFromForce tests how combinations of E, F, and v
% affect displacement u.

%% F vs. u when E=8000Pa with v=0.1 and 0.5
E=8000;
neu=0.1;
neu2=0.5;
forceType = 'groupForce';
nSample = 25;
stepForce = 20;
maxForce = 500;
nForce = round(maxForce/stepForce);
uxmat = zeros(nSample,nForce); % magnitudes of ux near the peak
uymat = zeros(nSample,nForce); % magnitudes of uy near the peak
uxmatIncomp = zeros(nSample,nForce); % magnitudes of ux near the peak
uymatIncomp = zeros(nSample,nForce); % magnitudes of uy near the peak
p=0;
parfor p=1:nForce
    f = stepForce*p;
    [ux,uy] = testDisplacementFromForce(forceType,f,E,neu);
    [uxIC,uyIC] = testDisplacementFromForce(forceType,f,E,neu2);
    [uys, uyID] = sort(uy(:));
    uymat(:,p) = abs(uys(1:nSample));
    uxmat(:,p) = abs(ux(uyID(1:nSample)));
    [uysIC, uyIDIC] = sort(uyIC(:));
    uymatIncomp(:,p) = abs(uysIC(1:nSample));
    uxmatIncomp(:,p) = abs(uxIC(uyIDIC(1:nSample)));
end
% %% plotting F vs. u
% figure, boxplot(uxmat,'plotstyle','compact','Color','g')
% hold on
% boxplot(uymat,'plotstyle','compact','Color','b')
% boxplot(uymatIncomp,'plotstyle','compact','Color','r')
% boxplot(uxmatIncomp,'plotstyle','compact','Color','c')
% plot(1:nForce',median(uymat),'b')
% plot(1:nForce',median(uymatIncomp),'r')
% plot(1:nForce',median(uxmat),'g')
% plot(1:nForce',median(uxmatIncomp),'c')
% ylim([0 max(uymat(:))*1.1])
% xlabel('Traction (Pa)')
% ylabel('Displacement u (pixel)')
% text(1,0.0155,'E=8kPa','FontSize',18)
% legend('u_y, \v=0.1','u_y, \v=0.5','u_x, \v=0.1','u_x, \v=0.5')
% % 'u_x, \v=0.1','u_y, \v=0.1','u_y, \v=0.5','u_x, \v=0.5',
%% plotting F vs. u
% figure, boxplot(uxmat,'notch','on','Color','g')
% hold on
% boxplot(uymat,'notch','on','Color','b')
% boxplot(uymatIncomp,'notch','on','Color','r')
% boxplot(uxmatIncomp,'notch','on','Color','c')
figure, hold on
errorbar(stepForce:stepForce:maxForce,mean(uymat),std(uymat),'-bs')
errorbar(stepForce:stepForce:maxForce,mean(uymatIncomp),std(uymatIncomp),'-rs')
errorbar(stepForce:stepForce:maxForce,mean(uxmat),std(uxmat),'-gs')
errorbar(stepForce:stepForce:maxForce,mean(uxmatIncomp),std(uxmatIncomp),'-cs')
ylim([0 max(uymat(:))*1.1])
xlabel('Traction (Pa)')
ylabel('Displacement u (pixel)')
text(stepForce,0.0155,'E = 8 kPa','FontSize',18)
legend('u_y,Neu=0.1','u_y,Neu=0.5','u_x,Neu=0.1','u_x,Neu=0.5')
xlim([0 maxForce])
%% u vs. E when F=100Pa
f=100;
nSample = 25;
stepE = 100;
maxE = 10000;
nE = round(maxE/stepE);
uxmatE = zeros(nSample,nE); % magnitudes of ux near the peak
uymatE = zeros(nSample,nE); % magnitudes of uy near the peak
uxmatIncompE = zeros(nSample,nE); % magnitudes of ux near the peak
uymatIncompE = zeros(nSample,nE); % magnitudes of uy near the peak
parfor p=1:nE
    E = stepE*p;
    [ux,uy] = testDisplacementFromForce(forceType,f,E,neu);
    [uxIC,uyIC] = testDisplacementFromForce(forceType,f,E,neu2);
    [uys, uyID] = sort(uy(:));
    uymatE(:,p) = abs(uys(1:nSample));
    uxmatE(:,p) = abs(ux(uyID(1:nSample)));
    [uysIC, uyIDIC] = sort(uyIC(:));
    uymatIncompE(:,p) = abs(uysIC(1:nSample));
    uxmatIncompE(:,p) = abs(uxIC(uyIDIC(1:nSample)));
end
%% plotting
% figure, boxplot(uxmatE,'notch','on','Color','g')
% hold on
% boxplot(uxmatIncompE,'notch','on','Color','c')
% boxplot(uymatE,'notch','on','Color','b')
% boxplot(uymatIncompE,'notch','on','Color','r')
figure, hold on
errorbar(stepE:stepE:maxE,mean(uymatE),std(uymatE),'-bs')
errorbar(stepE:stepE:maxE,mean(uymatIncompE),std(uymatIncompE),'-rs')
errorbar(stepE:stepE:maxE,mean(uxmatE),std(uxmatE),'-gs')
errorbar(stepE:stepE:maxE,mean(uxmatIncompE),std(uxmatIncompE),'-cs')
ylim([0 max(uymatE(:))*1.1])
xlabel('Elastic modulus (Pa)')
ylabel('Displacement u (pixel)')
text(6500,0.02,'F=100Pa','FontSize',18)
legend('uy, \Neu=0.1','uy, \Neu=0.5','ux, \Neu=0.1','ux, \Neu=0.5')
xlim([0 8000])
%% close up in low disp
ylim([0 0.025])
%% F vs. E
nSample = 25;
stepE = 200;
maxE = 5000;
nE = round(maxE/stepE);
stepForce = 100;
maxForce = 2000;
nForce = round(maxForce/stepForce);
uxmeanFE = zeros(nForce,nE); 
uymeanFE = zeros(nForce,nE); 
uxmeanFEIC = zeros(nForce,nE); 
uymeanFEIC = zeros(nForce,nE); 
h = waitbar(0,'Doing first iteration ...');
for p=1:nForce
    t= cputime;
    parfor q=1:nE
        E = stepE*q;
        f = stepForce*p;
        [ux,uy] = testDisplacementFromForce(forceType,f,E,neu);
        [uxIC,uyIC] = testDisplacementFromForce(forceType,f,E,neu2);
        [uys, uyID] = sort(uy(:));
        uymeanFE(p,q) = mean(abs(uys(1:nSample)));
        uxmeanFE(p,q) = mean(abs(ux(uyID(1:nSample))));
        [uysIC, uyIDIC] = sort(uyIC(:));
        uymeanFEIC(p,q) = mean(abs(uysIC(1:nSample)));
        uxmeanFEIC(p,q) = mean(abs(uxIC(uyIDIC(1:nSample))));
    end
    waitbar(p/nForce,h,[num2str((nForce-p)*(cputime-t)) 'sec left.'])
end

%% filling E=200:100:2400 vs. F=0:10:100
stepE2 = 100;
maxE2 = 2400;
nE2 = round(maxE2/stepE2);
stepForce2 = 10;
maxForce2 = 100;
nForce2 = round(maxForce2/stepForce2);
uxmeanFEfine = zeros(nForce2,nE2); 
uymeanFEfine = zeros(nForce2,nE2); 
uxmeanFEICfine = zeros(nForce2,nE2); 
uymeanFEICfine = zeros(nForce2,nE2); 
h2 = waitbar(0,'Doing first iteration ...');
for p=1:nForce2
    t= cputime;
    parfor q=1:nE2
        E = stepE2*q;
        f = stepForce2*p;
        [ux,uy] = testDisplacementFromForce(forceType,f,E,neu);
        [uxIC,uyIC] = testDisplacementFromForce(forceType,f,E,neu2);
        [uys, uyID] = sort(uy(:));
        uymeanFEfine(p,q) = mean(abs(uys(1:nSample)));
        uxmeanFEfine(p,q) = mean(abs(ux(uyID(1:nSample))));
        [uysIC, uyIDIC] = sort(uyIC(:));
        uymeanFEICfine(p,q) = mean(abs(uysIC(1:nSample)));
        uxmeanFEICfine(p,q) = mean(abs(uxIC(uyIDIC(1:nSample))));
    end
    waitbar(p/nForce,h2,[num2str((nForce-p)*(cputime-t)) 'sec left.'])
end
%% filling E=2400:1300:5000 vs. F=10:45:100
stepE3 = 1300;
maxE3 = 5000;
nE3 = 3;
stepForce3 = 45;
maxForce3 = 100;
nForce3 = 3;
uxmeanFEfine3 = zeros(nForce3,nE3); 
uymeanFEfine3 = zeros(nForce3,nE3); 
uxmeanFEICfine3 = zeros(nForce3,nE3); 
uymeanFEICfine3 = zeros(nForce3,nE3); 
h2 = waitbar(0,'Doing first iteration ...');
for p=1:nForce3
    t= cputime;
    parfor q=1:nE3
        E = stepE3*(q-1)+2400;
        f = stepForce3*(p-1)+10;
        [ux,uy] = testDisplacementFromForce(forceType,f,E,neu);
        [uxIC,uyIC] = testDisplacementFromForce(forceType,f,E,neu2);
        [uys, uyID] = sort(uy(:));
        uymeanFEfine3(p,q) = mean(abs(uys(1:nSample)));
        uxmeanFEfine3(p,q) = mean(abs(ux(uyID(1:nSample))));
        [uysIC, uyIDIC] = sort(uyIC(:));
        uymeanFEICfine3(p,q) = mean(abs(uysIC(1:nSample)));
        uxmeanFEICfine3(p,q) = mean(abs(uxIC(uyIDIC(1:nSample))));
    end
    waitbar(p/nForce3,h2,[num2str((nForce3-p)*(cputime-t)) 'sec left.'])
end
%% Visualizaing
f=stepForce:stepForce:maxForce;
E=stepE:stepE:maxE;
dataPath = '/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/testDisplacementFromForce';
% visualizeError(f,E,uymeanFE,dataPath,'pcolor_with_level1line',0.1)
figure,
[Egrid,Fgrid]=meshgrid(E,f);
Efine=stepE:stepE2:maxE;
ffine = stepForce:stepForce2:maxForce;
[Egridfine,Fgridfine]=meshgrid(Efine,ffine);
uymeanFEfine1 = interp2(Egrid,Fgrid,uymeanFE,Egridfine,Fgridfine);
pcolor(Efine,ffine,uymeanFEfine1)
shading interp

% pcolor(E,f,uymeanFE)
% shading interp
caxis([0 0.1])
colormap hot
colorbar('location','eastoutside')
hold on
contour(E,f,uymeanFE,[0.01 0.01],'LineWidth',1,'LineColor',[63/255,162/255,10/255]);
xlabel('Gel elastic modulus (Pa)')
ylabel('Traction magnitude (Pa)')
f2=stepForce2:stepForce2:maxForce2;
E2=stepE2:stepE2:maxE2;
hold on
pcolor(E2,f2,uymeanFEfine)
shading interp
contour(E2,f2,uymeanFEfine,[0.01 0.01],'LineWidth',1,'LineColor',[63/255,162/255,10/255]);
ylim([10 500])
E3=maxE2:stepE3:maxE;
f3 = stepForce2:stepForce3:maxForce2;
[Egrid3,Fgrid3]=meshgrid(E3,f3);
E3fine=maxE2:stepE2:maxE;
f3fine = stepForce2:stepForce2:maxForce2;
[Egrid3fine,Fgrid3fine]=meshgrid(E3fine,f3fine);
uymeanFE3 = interp2(Egrid3,Fgrid3,uymeanFEfine3,Egrid3fine,Fgrid3fine);
pcolor(E3fine,f3fine,uymeanFE3)
shading interp
colormap(flipud(colormap))


%% For incompressible gel
figure,
uymeanFEICfine1 = interp2(Egrid,Fgrid,uymeanFEIC,Egridfine,Fgridfine);
pcolor(Efine,ffine,uymeanFEICfine1)
shading interp
caxis([0 0.1])
colormap hot
colorbar('location','eastoutside')
hold on
contour(E,f,uymeanFEIC,[0.01 0.01],'LineWidth',1,'LineColor',[63/255,162/255,10/255]);
xlabel('Gel elastic modulus (Pa)')
ylabel('Traction magnitude (Pa)')
pcolor(E2,f2,uymeanFEICfine)
shading interp
contour(E2,f2,uymeanFEICfine,[0.01 0.01],'LineWidth',1,'LineColor',[63/255,162/255,10/255]);
ylim([10 500])

uymeanFEIC3 = interp2(Egrid3,Fgrid3,uymeanFEICfine3,Egrid3fine,Fgrid3fine);
pcolor(E3fine,f3fine,uymeanFEIC3)
shading interp
colormap(flipud(colormap))
%interpolation to blanked area
%% comparison between compressible vs. incompressible gel
figure,
contour(E,f,uymeanFEIC,[0.01 0.01],'LineWidth',1,'LineColor',[63/255,162/255,10/255]);
hold on
contour(E2,f2,uymeanFEICfine,[0.01 0.01],'LineWidth',1,'LineColor',[63/255,162/255,10/255]);
contour(E,f,uymeanFE,[0.01 0.01],'LineWidth',1,'LineColor','r');
contour(E2,f2,uymeanFEfine,[0.01 0.01],'LineWidth',1,'LineColor','r');
ylim([10 300])
