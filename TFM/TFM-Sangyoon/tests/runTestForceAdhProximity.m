% runTestForceAdhProximity runs testForceAdhProximity
%% Simulations - simple enough simulation - only d
nExp = 5;
f = 400; % Pa
r = 4; % pixel
expName = ['r' num2str(r)];
dmin = 2*r;
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
bead_x = cell(nExp,1);
bead_y = cell(nExp,1);
Av = cell(nExp,1);
cropInfoL2 = cell(nExp,nDist);
cropInfoL1 = cell(nExp,nDist);
for epm=1:nExp
    ii=0;
    for d=dmin:2:dmax % distance between adhesions
        ii=ii+1;
        if ii==1
            method = 'L2';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL2(epm,ii),nDetectedL2(epm,ii),fm1L2(epm,ii),fm2L2(epm,ii),fMapL2{epm,ii},cropInfoL2{epm,ii},bead_x{epm}, bead_y{epm}, Av{epm}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath);
            method = 'L1';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL1(epm,ii),nDetectedL1(epm,ii),fm1L1(epm,ii),fm2L1(epm,ii),fMapL1{epm,ii},cropInfoL1{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath,bead_x{epm}, bead_y{epm}, Av{epm});
        else
            method = 'L2';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL2(epm,ii),nDetectedL2(epm,ii),fm1L2(epm,ii),fm2L2(epm,ii),fMapL2{epm,ii},cropInfoL2{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath,bead_x{epm}, bead_y{epm}, Av{epm});
            method = 'L1';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL1(epm,ii),nDetectedL1(epm,ii),fm1L1(epm,ii),fm2L1(epm,ii),fMapL1{epm,ii},cropInfoL1{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath,bead_x{epm}, bead_y{epm}, Av{epm});
        end
    end
end
%% save and plotting
save(['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/' expName '/d_effect.mat'])
% dataPath = '/hms/scratch1/sh268/singleForceTesting/f_vs_d/forceDetec';

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
%% Figure making
% array of all reconstructed stresses
fh = figure;
tmax = 100;
midx = 51;
midy = 51;
for epm=1:nExp
    epmRev = nExp-epm+1;
    for ii=1:nDist % distance between adhesions
        axes('Position',[(ii-1)*(1/nDist), ((epmRev-1)*2+1)/(2*nExp), 1/nDist, 1/(2*nExp)])
        imshow(fMapL2{epm,ii},[0 tmax])
        colormap jet
        d = dmin+(ii-1)*2;
        hold on, plot(midx+d/2-cropInfoL2{epm,ii}(1)+1,midy-cropInfoL2{epm,ii}(2)+1,'ro')
        plot(midx-d/2-cropInfoL2{epm,ii}(1)+1,midy-cropInfoL2{epm,ii}(2)+1,'ro')
        axes('Position',[(ii-1)*(1/nDist), ((epmRev-1)*2)/(2*nExp), 1/nDist, 1/(2*nExp)])
        imshow(fMapL1{epm,ii},[0 tmax])        
        colormap jet
        hold on, plot(midx+d/2-cropInfoL1{epm,ii}(1)+1,midy-cropInfoL1{epm,ii}(2)+1,'ro')
        plot(midx-d/2-cropInfoL1{epm,ii}(1)+1,midy-cropInfoL1{epm,ii}(2)+1,'ro')
    end
end
%% save and plotting
save('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/d_effect.mat')

