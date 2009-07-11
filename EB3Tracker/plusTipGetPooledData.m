function [movGroup]=plusTipGetPooledData(mData,movGroup,prop2sample)


if iscell(prop2sample)
    nProps=max(size(prop2sample));
elseif isstr(prop2sample)
    nProps=1;
else
    error('plusTipGetPooledData: prop2sample must be string or cell array')
end

for iProp=1:nProps;
    if iscell(prop2sample)
        tempProp=prop2sample{iProp};
    else
        tempProp=prop2sample;
    end

    nGroups=length(movGroup);
    for iGroup=1:nGroups
        nMovs=length(movGroup(iGroup,1).movNums);
        pooledData=cell(nMovs,1);
        for iMov=1:nMovs

            % load projData for iMovie
            dataLoc=[mData{movGroup(iGroup,1).movNums(iMov,1),7}];
            load([dataLoc filesep 'meta' filesep 'projData']);

            % extract the track profile matrix
            tempProf=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
            switch tempProp
                case 'growthSpeeds'
                    % extract growth speeds
                    iMovData=tempProf(tempProf(:,5)==1,4);
                case 'shrinkSpeeds'
                    % extract shrinkage speeds
                    iMovData=tempProf(tempProf(:,5)==3,4);
                case 'gapLength'
                    % extract pause/shrinkage lengths (framees)
                    iMovData=[tempProf(tempProf(:,5)==2,6); tempProf(tempProf(:,5)==3,6)];
                case 'lifeTimes'
                    % extract lifetime distributions
                    lft=nan(projData.numTracks,1);
                    for iTrack=1:projData.numTracks
                        idx=find(tempProf(:,1)==iTrack);
                        if tempProf(idx(1),2)==1 | tempProf(idx(end),3)==projData.numFrames
                            %do nothing, track is at start or end of movie
                        else
                            lft(iTrack)=sum(tempProf(idx,6))+1;
                        end
                    end
                    lft(isnan(lft))=[];
                    iMovData=lft;

                otherwise
                    error('prop2sample not supported')
            end
            pooledData{iMov,1}=iMovData;
        end
        movGroup(iGroup,1).(tempProp)=cell2mat(pooledData);
    end
end