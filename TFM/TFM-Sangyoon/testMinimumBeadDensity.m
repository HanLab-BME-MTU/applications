% testMinimumBeadDensity calculates force detectability depending on bead
% density and white noise.
clear
%% Simulations - with pixel-wise tracking
numBeads = [250 500 750 1000 1500 2000 2500 3000 4000];
nNB = length(numBeads);
pixSize = 108; %nm/px
pixSizeUm = pixSize/1000;
d = 4; %pix (which means ~400 nm)
noiseLevels = 0:0.05:0.2;
nNL = length(noiseLevels);
nExp = 3;
E=2000; %Pa

d_err_Adh = zeros(nNB,nNL,nExp);
d_err_BG = zeros(nNB,nNL,nExp);
f_err_ADh = zeros(nNB,nNL,nExp);
f_err_BG = zeros(nNB,nNL,nExp);
dispDetec = zeros(nNB,nNL,nExp);
forceDetec = zeros(nNB,nNL,nExp);
pFR = zeros(nNB,nNL,nExp);
beadsOnAdh = zeros(nNB, nExp);
cL = 15;%[9 15 21]
adhArea = pi/4*(d*pixSize)^2*1e-18; %m2
f = (15e-12)/adhArea; %15 picoNewton/adhesion area, Pa
adhAreaUm2 = adhArea*1e12;
w = 2^8; h=w;
fieldArea = w*h; %pix
fieldAreaUm2 = fieldArea*pixSizeUm^2;
beadDensity = numBeads/fieldAreaUm2; % #/um2
B2Bspacing = 1./beadDensity.^2;
% kk=0;
%% simulation for numBeads and noiseLevel
pathBasisClassTbl='/storage/network/TFM_Development/TFM2D/TFM_basisClass/basisClass2kPa9pix.mat';
for epm=1:nExp
    p=0;
    ii=0;
    for nB=numBeads
        ii=ii+1;
        jj=0;
        for curWN = noiseLevels
            jj=jj+1;
            p=p+1;
            dataPath=['/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/f' num2str(round(f)) ...
                'd' num2str(d) 'nB' num2str(nB) 'noise' num2str(curWN) 'exp' num2str(epm)];

            if jj==1
                [d_err_Adh(ii,jj,epm),d_err_BG(ii,jj,epm),dispDetec(ii,jj,epm),...
                    f_err_ADh(ii,jj,epm),f_err_BG(ii,jj,epm),pFR(ii,jj,epm),forceDetec(ii,jj,epm),...
                    beadsOnAdh(ii,epm),bead_x{ii,epm}, bead_y{ii,epm}, Av{ii,epm}]= ...
                    testSingleForce(f,d,cL,dataPath,[],[],[],'QR','whiteNoise',curWN,'nPoints',nB,...
                    'pathBasisClassTbl',pathBasisClassTbl,'E',E);
            else
                [d_err_Adh(ii,jj,epm),d_err_BG(ii,jj,epm),dispDetec(ii,jj,epm),...
                    f_err_ADh(ii,jj,epm),f_err_BG(ii,jj,epm),pFR(ii,jj,epm),forceDetec(ii,jj,epm),...
                    beadsOnAdh(ii,epm)]= ...
                    testSingleForce(f,d,cL,dataPath,bead_x{ii,epm}, bead_y{ii,epm}, Av{ii,epm},...
                    'QR','whiteNoise',curWN, 'pathBasisClassTbl',pathBasisClassTbl,'E',E);
            end
        end
    end
end
%% plotting
myCs=distinguishable_colors(nNL);
figure, hold on
for jj=1:nNL
    curCol = myCs(jj,:);
    brighterCol = curCol.*1.5;
    brighterCol(brighterCol>1) = 1;
    plot(beadDensity,squeeze(forceDetec(:,jj,:)),'Color',brighterCol);
    h(jj) = plot(beadDensity,mean(forceDetec(:,jj,:),3),'LineWidth',2,'Color',myCs(jj,:));
end
legend(h,'0 %','5 %','10%','15%','20%')
savefig('/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/beadDensityVsForceDetec.fig')
%% save
dataPath='/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/allData.mat';
save(dataPath)
