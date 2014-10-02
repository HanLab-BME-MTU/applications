% runTestForceAdhProximity runs testForceAdhProximity
%% Simulations - simple enough simulation - only d
nExp = 30;
f = 400; % Pa
r = 2; % radius in pixel
expName = ['dia' num2str(2*r)];
dmin = 2*r; % center-to-center distance between the two adhesions
dmax = 30; % pixel
nDist = (dmax-dmin)/2+1;
detectedL2 = zeros(nExp,nDist);
detectedL1 = zeros(nExp,nDist);
nDetectedL2 = zeros(nExp,nDist);
nDetectedL1 = zeros(nExp,nDist);
fm1L2 = zeros(nExp,nDist);
fm2L2 = zeros(nExp,nDist);
fm1L1 = zeros(nExp,nDist);
fm2L1 = zeros(nExp,nDist);
fMapL2 = cell(nExp,nDist);
fMapL1 = cell(nExp,nDist);
oMapL2 = cell(nExp,nDist);
oMapL1 = cell(nExp,nDist);
bead_x = cell(nExp,1);
bead_y = cell(nExp,1);
Av = cell(nExp,1);
cropInfoL2 = cell(nExp,nDist);
cropInfoL1 = cell(nExp,nDist);
%% Simulations - simple enough simulation - only d
for epm=1:nExp
    ii=0;
    for d=dmin:2:dmax % distance between adhesions
        ii=ii+1;
        if ii==1
            method = 'L2';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL2(epm,ii),nDetectedL2(epm,ii),fm1L2(epm,ii),fm2L2(epm,ii),fMapL2{epm,ii},cropInfoL2{epm,ii},bead_x{epm}, bead_y{epm}, Av{epm}, oMapL2{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath);
            method = 'L1';
            oldDataPath = dataPath;
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL1(epm,ii),nDetectedL1(epm,ii),fm1L1(epm,ii),fm2L1(epm,ii),fMapL1{epm,ii},cropInfoL1{epm,ii},~,~,~, oMapL2{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath,oldDataPath);
        else
            method = 'L2';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL2(epm,ii),nDetectedL2(epm,ii),fm1L2(epm,ii),fm2L2(epm,ii),fMapL2{epm,ii},cropInfoL2{epm,ii},~,~,~, oMapL2{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath,[],bead_x{epm}, bead_y{epm}, Av{epm});
            oldDataPath = dataPath;
            method = 'L1';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL1(epm,ii),nDetectedL1(epm,ii),fm1L1(epm,ii),fm2L1(epm,ii),fMapL1{epm,ii},cropInfoL1{epm,ii},~,~,~, oMapL2{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath,oldDataPath);
        end
    end
end
%% For only left adhesion
fm1L2left = zeros(nExp,nDist);
fm2L2left = zeros(nExp,nDist);
fm1L1left = zeros(nExp,nDist);
fm2L1left = zeros(nExp,nDist);
fMapL2left = cell(nExp,nDist);
fMapL1left = cell(nExp,nDist);
cropInfoL2left = cell(nExp,nDist);
cropInfoL1left = cell(nExp,nDist);
oMapL2left = cell(nExp,nDist);
oMapL1left = cell(nExp,nDist);
for epm=1:nExp
    ii=0;
    for d=dmin:2:dmax % distance between adhesions
        ii=ii+1;
        method = 'L2';
        dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method 'One'];
        [~,~,fm1L2left(epm,ii),fm2L2left(epm,ii),fMapL2left{epm,ii},cropInfoL2left{epm,ii},~,~,~,oMapL2left{epm,ii}] = ...
            testForceAdhProximity(d,f,0,r,r,method,dataPath,[],bead_x{epm}, bead_y{epm}, Av{epm});
        oldDataPath = dataPath;
        method = 'L1';
        dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method 'One'];
        [~,~,fm1L1left(epm,ii),fm2L1left(epm,ii),fMapL1left{epm,ii},cropInfoL1left{epm,ii},~,~,~,oMapL1left{epm,ii}] = ...
            testForceAdhProximity(d,f,0,r,r,method,dataPath,oldDataPath);
    end
end

%% Plotting
figure, plot(dmin:2:dmax,mean(detectedL2))
hold all
ylim([0 1.2])
plot(dmin:2:dmax,mean(detectedL1))
title('detected')
%% The number of detected force loc max
figure, plot(dmin:2:dmax,mean(nDetectedL2))
hold all
ylim([0 1.2])
plot(dmin:2:dmax,mean(nDetectedL1))
title('the number of detected force local maxima')
%% force mag
figure, plot(dmin:2:dmax,mean([fm1L2; fm2L2]))
hold all
plot(dmin:2:dmax,mean([fm1L1; fm2L1]))
title('force magnitude reconstructed')
% %% recalculation of detected
% for epm=1:nExp
%     ii=0;
%     for d=dmin:2:dmax % distance between adhesions
%         ii=ii+1;
%         x1 = 51-d/2;
%         y1=51;
%         y2=51;
%         x2 = 51+d/2;
%         [sizeL2y, sizeL2x] = size(fMapL2{epm,ii});
%         [sizeL1y, sizeL1x] = size(fMapL1{epm,ii});
%         cropInfoL2{epm,ii} = [(100-sizeL2x)/2 (100-sizeL2y)/2];
%         cropInfoL1{epm,ii} = [(100-sizeL1x)/2 (100-sizeL1y)/2];
%         pstructL2 = pointSourceDetection(fMapL2{epm,ii},1.1);
%         fMaximaL2 = [pstructL2.x'+cropInfoL2{epm,ii}(1)-1 pstructL2.y'+cropInfoL2{epm,ii}(2)-1];
%         pstructL1 = pointSourceDetection(fMapL1{epm,ii},1.1);
%         fMaximaL1 = [pstructL1.x'+cropInfoL1{epm,ii}(1)-1 pstructL1.y'+cropInfoL1{epm,ii}(2)-1];
%         if ~isempty(pstructL2)
%             [~,dist1L2] = KDTreeClosestPoint(fMaximaL2, [51-d/2 51]);
%             [~,dist2L2] = KDTreeClosestPoint(fMaximaL2, [51+d/2 51]);
%             if dist1L2<7 && dist2L2<7
%                 detectedL2(epm,ii) = 1;
%             elseif dist1L2<7 || dist2L2<7
%                 detectedL2(epm,ii) = 0.5;
%             else
%                 detectedL2(epm,ii) = 0;
%             end
%         else
%             detectedL2(epm,ii) = 0;
%         end
%         if ~isempty(pstructL1)
%             [~,dist1L1] = KDTreeClosestPoint(fMaximaL1, [51-d/2 51]);
%             [~,dist2L1] = KDTreeClosestPoint(fMaximaL1, [51+d/2 51]);
%             if dist1L1<7 && dist2L1<7
%                 detectedL1(epm,ii) = 1;
%             elseif dist1L1<7 || dist2L1<7
%                 detectedL1(epm,ii) = 0.5;
%             else
%                 detectedL1(epm,ii) = 0;
%             end
%         else
%             detectedL1(epm,ii) = 0;
%         end
%     end
% end
%% Figure making - two adhesions
% array of all reconstructed stresses
fh = figure;
tmax = 100;
midx = 51;
midy = 51;
colormap jet
skippingN = 6;
actualNexp = floor(nExp/skippingN);
for ii=1:nDist
    axes('Position',[(ii-1)*(1/nDist), (2*actualNexp)/(2*actualNexp+1), 1/nDist, 1/(2*actualNexp+1)])
    imshow(oMapL2{1,ii},[0 400])
end
for epm=skippingN:skippingN:nExp%1:nExp
    epmRev = actualNexp-floor(epm/skippingN)+1;
    for ii=1:nDist % distance between adhesions
        axes('Position',[(ii-1)*(1/nDist), ((epmRev-1)*2+1)/(2*actualNexp+1), 1/nDist, 1/(2*actualNexp+1)])
        imshow(fMapL2{epm,ii},[0 tmax])
        colormap jet
        d = dmin+(ii-1)*2;
        hold on, plot(midx+d/2-cropInfoL2{epm,ii}(1)+1,midy-cropInfoL2{epm,ii}(2)+1,'ro')
        plot(midx-d/2-cropInfoL2{epm,ii}(1)+1,midy-cropInfoL2{epm,ii}(2)+1,'ro')
        axes('Position',[(ii-1)*(1/nDist), ((epmRev-1)*2)/(2*actualNexp+1), 1/nDist, 1/(2*actualNexp+1)])
        imshow(fMapL1{epm,ii},[0 tmax])        
        colormap jet
        hold on, plot(midx+d/2-cropInfoL1{epm,ii}(1)+1,midy-cropInfoL1{epm,ii}(2)+1,'ro')
        plot(midx-d/2-cropInfoL1{epm,ii}(1)+1,midy-cropInfoL1{epm,ii}(2)+1,'ro')
    end
end
hgsave(fh,['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/twoAdhForceMap.fig'],'-v7.3')
print(fh,'-depsc2', '-r150', ['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/twoAdhForceMap.eps']);
%% Figure making - one adhesion
% array of all reconstructed stresses
fh2 = figure;

for ii=1:nDist
    axes('Position',[(ii-1)*(1/nDist), (2*actualNexp)/(2*actualNexp+1), 1/nDist, 1/(2*actualNexp+1)])
    imshow(oMapL2left{1,ii},[0 400])
end

for epm=skippingN:skippingN:nExp%1:nExp
    epmRev = actualNexp-floor(epm/skippingN)+1;%epmRev = nExp-epm+1;
    for ii=1:nDist % distance between adhesions
        axes('Position',[(ii-1)*(1/nDist), ((epmRev-1)*2+1)/(2*actualNexp+1), 1/nDist, 1/(2*actualNexp+1)])
        imshow(fMapL2left{epm,ii},[0 tmax])
        colormap jet
        d = dmin+(ii-1)*2;
        hold on, plot(midx+d/2-cropInfoL2left{epm,ii}(1)+1,midy-cropInfoL2left{epm,ii}(2)+1,'yo')
        plot(midx-d/2-cropInfoL2left{epm,ii}(1)+1,midy-cropInfoL2left{epm,ii}(2)+1,'ro')
        axes('Position',[(ii-1)*(1/nDist), ((epmRev-1)*2)/(2*actualNexp+1), 1/nDist, 1/(2*actualNexp+1)])
        imshow(fMapL1left{epm,ii},[0 tmax])        
        colormap jet
        hold on, plot(midx+d/2-cropInfoL1left{epm,ii}(1)+1,midy-cropInfoL1left{epm,ii}(2)+1,'yo')
        plot(midx-d/2-cropInfoL1left{epm,ii}(1)+1,midy-cropInfoL1left{epm,ii}(2)+1,'ro')
    end
end
hgsave(fh2,['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/oneAdhForceMap.fig'],'-v7.3')
print(fh2,'-depsc2', '-r150', ['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/oneAdhForceMap.eps']);

%% Ratio plot R = fm2/fm2left
ratioL2 = fm2L2left./fm2L2;
ratioL1 = fm2L1left./fm2L1;
fh3=figure; 
plot(dmin:2:dmax,ratioL2)
hold on
plot(dmin:2:dmax,mean(ratioL2))
plot(dmin:2:dmax,mean(ratioL1),'g')
hgsave(fh3,['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/ratio.fig'],'-v7.3')
print(fh3,'-depsc2', '-r150', ['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/ratio.eps']);
%% Plotting fm2 and fm2left -L2
fh4=figure; 
plot(dmin:2:dmax,fm2L2left)
hold on
plot(dmin:2:dmax,fm2L2,'k')
plot(dmin:2:dmax,mean(fm2L2left),'r', 'LineWidth',3)
plot(dmin:2:dmax,mean(fm2L2),'k','LineWidth',3)

hgsave(fh4,['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/fm2L2.fig'],'-v7.3')
print(fh4,'-depsc2', '-r150', ['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/fm2L2.eps']);
%% Plotting fm2 and fm2left -L1
fh5=figure; 
plot(dmin:2:dmax,fm2L1left , 'r')
hold on
plot(dmin:2:dmax,fm2L1,'k')
plot(dmin:2:dmax,mean(fm2L1left),'r', 'LineWidth',3)
plot(dmin:2:dmax,mean(fm2L1),'k','LineWidth',3)

hgsave(fh5,['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/fm2L1.fig'],'-v7.3')
print(fh5,'-depsc2', '-r150', ['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/fm2L1.eps']);
%% Plotting fm1 and fm1left - negative control - L2
fh6=figure; 
plot(dmin:2:dmax,fm1L2left,'r')
hold on
plot(dmin:2:dmax,fm1L2,'k')
plot(dmin:2:dmax,mean(fm1L2left),'r', 'LineWidth',3)
plot(dmin:2:dmax,mean(fm1L2),'k','LineWidth',3)

hgsave(fh6,['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/fm1L2.fig'],'-v7.3')
print(fh6,'-depsc2', '-r150', ['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/fm1L2.eps']);
%% Plotting fm1 and fm1left - negative control - L1
fh6=figure; 
plot(dmin:2:dmax,fm1L1left,'r')
hold on
plot(dmin:2:dmax,fm1L1,'k')
plot(dmin:2:dmax,mean(fm1L1left),'r', 'LineWidth',3)
plot(dmin:2:dmax,mean(fm1L1),'k','LineWidth',3)

hgsave(fh6,['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/fm1L1.fig'],'-v7.3')
print(fh6,'-depsc2', '-r150', ['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/fm1L1.eps']);
%% save 
fm2L2leftMean = mean(fm2L2left);
fm2L2Mean=mean(fm2L2);
fm2L1leftMean = mean(fm2L1left);
fm2L1Mean = mean(fm2L1);
fm1L2leftMean = mean(fm1L2left);
fm1L2Mean = mean(fm1L2);
fm1L1leftMean = mean(fm1L1left);
fm1L1Mean = mean(fm1L1);
save(['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/d_effect.mat'])
% dataPath = '/hms/scratch1/sh268/singleForceTesting/f_vs_d/forceDetec';
