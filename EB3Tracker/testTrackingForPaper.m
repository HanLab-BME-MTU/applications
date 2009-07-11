function [allData]=testTrackingForPaper


paramRange=2:2:50;
nIter=length(paramRange);
s1 = 2;
strg1 = sprintf('%%.%dd',s1);
for i=1:nIter
    timeWindow=paramRange(i);
    indxStr1 = sprintf(strg1,timeWindow);
    projList(i,1).imDir='X:\forKathryn\forPaperTest\ETH-RPE01-control-movie4\images';
    projList(i,1).anDir=['X:\forKathryn\forPaperTest\ETH-RPE01-control-movie4\roi_1_timeWin' filesep indxStr1 '_meta'];

    movGroup(i,1).movNums=i;
    movGroup(i,1).label=['ETH ' indxStr1];
end
save('projListETHtime','projList')
save('movGroupETHtime','movGroup')



allData=[];

testType=2;

if testType==1 % test time window parameter
    projList(1,1).imDir='X:\forKathryn\forPaperTest\ETH-RPE01-control-movie4\images';
    projList(1,1).anDir='X:\forKathryn\forPaperTest\ETH-RPE01-control-movie4\roi_1_timeWin';

    projList(2,1).imDir='X:\forKathryn\forPaperTest\NIH-HUVEC-control-movie4\images';
    projList(2,1).anDir='X:\forKathryn\forPaperTest\NIH-HUVEC-control-movie4\roi_1_timeWin';

    for iProj=1:2

        if iProj==1 %ETH

            %timeWindow=25;
            minTrackLen=3;

            minRadius=3;
            maxRadius=5;

            maxFAngle=45;
            maxShrinkFactor=1.5;
            d1Max=1.5;

            secPerFrame=0.8;
            pixSizeNm=110;
        else        %NIH

            %timeWindow=25;
            minTrackLen=3;

            minRadius=5;
            maxRadius=20;

            maxFAngle=45;
            maxShrinkFactor=1.5;
            d1Max=1.5;


            secPerFrame=2.0;
            pixSizeNm=105;
        end

        paramRange=2:2:50;
        nIter=length(paramRange);
        s1 = 2;
        strg1 = sprintf('%%.%dd',s1);
        data=cell(nIter,11);
        for iParam=1:nIter
            timeWindow=paramRange(iParam);

            plusTipTracker(projList(iProj),timeWindow,minTrackLen,minRadius,maxRadius,maxFAngle,maxShrinkFactor,d1Max);
            [projData]=plusTipPostTracking(projList(iProj),secPerFrame,pixSizeNm);

            indxStr1 = sprintf(strg1,timeWindow);
            movefile([projData.anDir filesep 'track'],[projData.anDir filesep indxStr1 '_track'])
            movefile([projData.anDir filesep 'meta'],[projData.anDir filesep indxStr1 '_meta'])

            temp=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;


            data{iParam,1}=timeWindow;
            data{iParam,2}=timeWindow*secPerFrame;
            data{iParam,3}=projData.numTracks;

            % calculate lifetimes of all tracks starting after first frame and
            % ending before last frame
            lft=nan(projData.numTracks,1);
            for iTrack=1:projData.numTracks
                idx=find(temp(:,1)==iTrack);
                if temp(idx(1),2)==1 | temp(idx(end),3)==projData.numFrames
                    %do nothing, track is at start or end of movie
                else
                    lft(iTrack)=sum(temp(idx,6))+1;
                end
            end
            data{iParam,4}=nanmean(lft);

            % get nPauses, avgPauseLength, and avgPauseDisplacement
            idx=find(temp(:,5)==2);
            gapLength=secPerFrame*(temp(idx,6)+1);
            dispMic=(pixSizeNm*temp(idx,7))/1000;

            data{iParam,5}=length(idx);
            data{iParam,6}=mean(gapLength);
            data{iParam,7}=mean(dispMic);

            % get nCat, avgShrinkLength, and avgShrinkDisplacement
            idx=find(temp(:,5)==3);
            gapLength=secPerFrame*(temp(idx,6)+1);
            dispMic=(pixSizeNm*temp(idx,7))/1000;

            data{iParam,8}=length(idx);
            data{iParam,9}=mean(gapLength);
            data{iParam,10}=abs(mean(dispMic));

            data{iParam,11}=projData.typeStats.type1_Ppause;
            data{iParam,12}=projData.typeStats.type1_Pcat;
        end

        allData(iProj,1).timeWindow=data;

    end
end




if testType==2 % test search radius range parameter
    projList(1,1).imDir='X:\forKathryn\forPaperTest\ETH-RPE01-control-movie4\images';
    projList(1,1).anDir='X:\forKathryn\forPaperTest\ETH-RPE01-control-movie4\roi_2_searchRad';

    projList(2,1).imDir='X:\forKathryn\forPaperTest\NIH-HUVEC-control-movie4\images';
    projList(2,1).anDir='X:\forKathryn\forPaperTest\NIH-HUVEC-control-movie4\roi_1_searchRad';

    for iProj=1:2

        if iProj==1 %ETH

            timeWindow=16;
            minTrackLen=3;

            minRadius=3.*ones(1,10);
            maxRadius=3:12;

            maxFAngle=45;
            maxShrinkFactor=1.5;
            d1Max=1.5;

            secPerFrame=0.8;
            pixSizeNm=110;
        else        %NIH

            timeWindow=16;
            minTrackLen=3;

            minRadius=8.*ones(1,10);
            maxRadius=8:2:26;

            maxFAngle=45;
            maxShrinkFactor=1.5;
            d1Max=1.5;


            secPerFrame=2.0;
            pixSizeNm=105;
        end

        nIter=length(minRadius);
        s1 = 2;
        strg1 = sprintf('%%.%dd',s1);
        data=cell(nIter,11);
        for iParam=1:nIter
            
            plusTipTracker(projList(iProj),timeWindow,minTrackLen,minRadius(iParam),maxRadius(iParam),maxFAngle,maxShrinkFactor,d1Max);
            [projData]=plusTipPostTracking(projList(iProj),secPerFrame,pixSizeNm);

            indxStr1 = sprintf(strg1,maxRadius(iParam));
            movefile([projData.anDir filesep 'track'],[projData.anDir filesep indxStr1 '_track'])
            movefile([projData.anDir filesep 'meta'],[projData.anDir filesep indxStr1 '_meta'])

            temp=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;

            close(1:gcf-1)
            saveas(gcf,[projData.anDir filesep indxStr1 '_meta' filesep 'searchRadMax_' indxStr1 '.fig'])
            saveas(gcf,[projData.anDir filesep 'searchRadMax_' indxStr1 '.png'])

            data{iParam,1}=minRadius(iParam);
            data{iParam,2}=maxRadius(iParam);
            data{iParam,3}=projData.numTracks;


        end

        allData(iProj,1).searchRad=data;

    end
end




