% scriptCheckFateSpatialIndependce
% load('movieData.mat')
load('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130429_cell11_100f/ROI/Colocalization/Region1_gc5/data/tracksNA.mat')
MD=MovieData.load('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130429_cell11_100f/ROI/movieData.mat');
% Do it for 130511 Cell2 before FAKi too
load('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130511 Cell2 beforeFAKi/ROItop/Colocalization/tracksSeparation/data/tracksNA.mat')
MD=MovieData.load('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130511 Cell2 beforeFAKi/ROItop/movieData.mat');
%% separateMatureAdhesionTracks
outputPath = [MD.getPath filesep 'Colocalization' filesep 'tracksSeparation'];
tracksNA = separateMatureAdhesionTracks(tracksNA, outputPath);
%% checkAdhesionFateSpatialIndependence
pixSize = MD.pixelSize_; % 72 nm/pixel
h1 = figure;
hold all
p=0;
for r_um = [2,4,6,8,16]
    p=p+1;
    r=r_um*1000/pixSize;
    [popStat(p),failNeiStat(p),matureNeiStat(p)] = checkAdhesionFateSpatialIndependence(tracksNA, r, outputPath);
    errorbar([nanmean(popStat(p).mRatio),nanmean(failNeiStat(p).mRatio),nanmean(matureNeiStat(p).mRatio)],...
                                [nanstd(popStat(p).mRatio)/sqrt(length(popStat(p).mRatio)),nanstd(failNeiStat(p).mRatio)/sqrt(length(failNeiStat(p).mRatio)),...
                                nanstd(matureNeiStat(p).mRatio)/sqrt(length(matureNeiStat(p).mRatio))])
end
% ylim([0 0.60])
r_um = [2,4,6,8,16];
ylabel('ratio of maturing NAs over all NAs at the emerging time')
legend(num2str(r_um(1)), num2str(r_um(2)),num2str(r_um(3)),num2str(r_um(4)),num2str(r_um(5)))
%% Statistical test
[h,p]=ttest2(failNeiStat(2).mRatio,matureNeiStat(2).mRatio)
[h,p]=ttest2(failNeiStat(3).mRatio,matureNeiStat(3).mRatio)
[h,p]=ttest2(failNeiStat(4).mRatio,matureNeiStat(4).mRatio)
[h,p]=ttest2(failNeiStat(5).mRatio,matureNeiStat(5).mRatio)

%% N
h2 = figure;
hold all
for r_um = [1,2,4,6,8,16,32]
    r=r_um*1000/pixSize;
    [popStat,failNeiStat,matureNeiStat] = checkAdhesionFateSpatialIndependence(tracksNA, r, outputPath);
    errorbar([mean(popStat.nFail+popStat.nMature),mean(failNeiStat.nMature+failNeiStat.nFail),mean(matureNeiStat.nFail+matureNeiStat.nMature)],...
                                [std(popStat.nFail+popStat.nMature)/sqrt(length(popStat.mRatio)),std(failNeiStat.nMature+failNeiStat.nFail)/sqrt(length(failNeiStat.mRatio)),...
                                std(matureNeiStat.nFail+matureNeiStat.nMature)/sqrt(length(matureNeiStat.mRatio))])
end
r_um = [1,2,4,6,8,16,32];
ylabel('the number')
% legend(num2str(r_um(1)), num2str(r_um(2)),num2str(r_um(3)),num2str(r_um(4)),num2str(r_um(5)))
legend(num2str(r_um(1)), num2str(r_um(2)),num2str(r_um(3)),num2str(r_um(4)),num2str(r_um(5)),num2str(r_um(6)),num2str(r_um(7)))

%% Showing maturing vs. failing adhesion tracks on cell image in different colors
% showing cell channel
matureID = arrayfun(@(x) (x.maturing==1),tracksNA);
failID = arrayfun(@(x) (x.maturing==0),tracksNA);
FAID = arrayfun(@(x) (x.maturing==3),tracksNA);
longNAID = arrayfun(@(x) (x.maturing==2),tracksNA);

matureIDidx = find(matureID);
failIDidx = find(failID);
FAIDidx = find(FAID);
longNAIDidx = find(longNAID);

paxChannel = MD.getChannel(2);
% For focal adhesions (maturing= 3) or nonmaturing long NA(maturing= 2) 
%% Showing itself
load('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130429_cell11_100f/ROI/Colocalization/trackSeparation/data/cropInfo.mat')
cropInfo = cropPosition;
f50 = 100;
figure, imshow(paxChannel.loadImage(f50),[])
hold on
for ii=matureIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        plot(tracksNA(ii).xCoord(tracksNA(ii).startingFrame:f50),tracksNA(ii).yCoord(tracksNA(ii).startingFrame:f50), 'g')
    end
end
for ii=failIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'r')
    end
end
for ii=FAIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        adhBidx = find(cellfun(@(x) ~isempty(x), tracksNA(ii).adhBoundary));
        if adhBidx(end) < f50 && adhBidx(end)>f50-10
            plot(tracksNA(ii).adhBoundary{adhBidx(end)}(:,2)+cropInfo(1), tracksNA(ii).adhBoundary{adhBidx(end)}(:,1)+cropInfo(2), 'k')
        elseif adhBidx(end) > f50
            % find the index that's closest to f50
            idxClosest = find(adhBidx<=f50,1,'last');
            plot(tracksNA(ii).adhBoundary{adhBidx(idxClosest)}(:,2)+cropInfo(1), tracksNA(ii).adhBoundary{adhBidx(idxClosest)}(:,1)+cropInfo(2), 'k')
        end
        plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'k')
    end
end
for ii=longNAIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'y')
    end
end

%% overlay with TFM
iiformat = ['%.' '3' 'd'];
forceImg = imread(['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130429_cell11_100f/ROI/Colocalization/trackSeparation/fMap/force' num2str(f50,iiformat) 'max1600.tif']);
figure, imshow(double(forceImg)*3500.0/2.0^15,[0 1200]), colormap jet
hold on
for ii=matureIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        plot(tracksNA(ii).xCoord(tracksNA(ii).startingFrame:f50)-cropInfo(1),tracksNA(ii).yCoord(tracksNA(ii).startingFrame:f50)-cropInfo(2), 'g')
    end
end
for ii=failIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        plot(tracksNA(ii).xCoord(1:f50)-cropInfo(1),tracksNA(ii).yCoord(1:f50)-cropInfo(2), 'r')
    end
end
for ii=FAIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        adhBidx = find(cellfun(@(x) ~isempty(x), tracksNA(ii).adhBoundary));
        if adhBidx(end) < f50 && adhBidx(end)>f50-10
            plot(tracksNA(ii).adhBoundary{adhBidx(end)}(:,2), tracksNA(ii).adhBoundary{adhBidx(end)}(:,1), 'k')
        elseif adhBidx(end) > f50
            % find the index that's closest to f50
            idxClosest = find(adhBidx<=f50,1,'last');
            plot(tracksNA(ii).adhBoundary{adhBidx(idxClosest)}(:,2), tracksNA(ii).adhBoundary{adhBidx(idxClosest)}(:,1), 'k')
        end
        plot(tracksNA(ii).xCoord(1:f50)-cropInfo(1),tracksNA(ii).yCoord(1:f50)-cropInfo(1), 'k')
    end
end
for ii=longNAIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        plot(tracksNA(ii).xCoord(1:f50)-cropInfo(1),tracksNA(ii).yCoord(1:f50)-cropInfo(1), 'y')
    end
end




