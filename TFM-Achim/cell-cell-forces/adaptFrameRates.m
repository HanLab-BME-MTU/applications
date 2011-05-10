function corrSets=adaptFrameRates(corrSets,tBtwFrms,aveType)
display('Adapted frame rates')
%**************************************************************************
% 2) Bring all time courses to the same frame rate                        *
%**************************************************************************

% check if the given corrSets are cells or edges:
if isfield(corrSets,'resF')
    isCells=1;
else
    isCells=0;
end

% These are either experiments with unmeasured dt or with super fast frame
% rates:
dtList = vertcat(corrSets.dt_mean);
badSets=dtList<=2;
% disregard those:
display(['# Disregarded experiments:= ',num2str(sum(badSets)),' out of:= ',num2str(length(corrSets))]);
corrSets(badSets)=[];
dtList(badSets)=[];
% find minimal dt:
dt_min = min(dtList);
% Warn the user if dt_min<tBtwFrms:
if dt_min<tBtwFrms
    display('The fastest frame rate is faster then the set time between frames!')
end

% Find the max time span:
t_max   = max(vertcat(corrSets.t));

% general time vector to all experiments will be mapped:
tvec=(0:tBtwFrms:(t_max+tBtwFrms))';

for idx=1:length(corrSets)
    % rounded time vector of current experiment:
    idsTPts=round(corrSets(idx).t./tBtwFrms)+1;
    intTPts=round(corrSets(idx).t./tBtwFrms)*tBtwFrms;
    % disregard frames that are bad anyways:
    goodVals=~isnan(idsTPts);
    idsTPts=idsTPts(goodVals);
    intTPts=intTPts(goodVals);
    
    %!!!% ToDo: check if there are multiple entries in idsTPts!
    
    if sum(tvec(idsTPts)-intTPts)~=0
        error('Check sum should be zero. Time points are not calculated correctly')
    end
    
    % build the new time vector:
    newTvec=NaN*zeros(size(tvec));
    newTvec(idsTPts)=intTPts;
    
    if nansum(tvec-newTvec)~=0
        error('Check sum should be zero. Timepts are not calculated correctly')
    end
    
    if isCells
        % first adapt the residual force:
        % do average if necessary:
        if strcmpi(aveType,'mean') || strcmpi(aveType,'nanmean')
            resF=aveInSameBin(corrSets(idx).resF,idsTPts,aveType);
        else
            resF=corrSets(idx).resF;
        end
        corrSetsNew(idx).resF = NaN*zeros(size(tvec,1),2);
        corrSetsNew(idx).resF(idsTPts,:) = resF(goodVals,:);
        
        % run through all edges and set the values at the right position:
        for edgeId=1:length(corrSets(idx).edge)
            if strcmpi(aveType,'mean') || strcmpi(aveType,'nanmean')
                fcorr=aveInSameBin(corrSets(idx).edge(edgeId).fcorr,idsTPts,aveType);
            else
                fcorr=corrSets(idx).edge(edgeId).fcorr;
            end
            corrSetsNew(idx).edge(edgeId).fcorr = NaN*zeros(size(tvec,1),2);
            corrSetsNew(idx).edge(edgeId).fcorr(idsTPts,:) = fcorr(goodVals,:);
        end
    else
        if strcmpi(aveType,'mean') || strcmpi(aveType,'nanmean')
            fcorr=aveInSameBin(corrSets(idx).fcorr,idsTPts,aveType);
            Icorr=aveInSameBin(corrSets(idx).Icorr,idsTPts,aveType);
        else
            fcorr=corrSets(idx).fcorr;
            Icorr=corrSets(idx).Icorr;
        end
        corrSetsNew(idx).fcorr = NaN*zeros(size(tvec,1),1);
        corrSetsNew(idx).Icorr = NaN*zeros(size(tvec,1),1);
        
        corrSetsNew(idx).fcorr(idsTPts,1) = fcorr(goodVals,:);
        corrSetsNew(idx).Icorr(idsTPts,1) = Icorr(goodVals,:);
    end
        
    % fill in the new time values:
    corrSetsNew(idx).t       = newTvec;
    corrSetsNew(idx).dt_mean = tBtwFrms;
    corrSetsNew(idx).dt_std  = std(corrSets(idx).t(goodVals)-intTPts);
end
% swith back to the old variable name:
clear corrSets
corrSets=corrSetsNew;