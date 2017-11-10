function corrSets=pullCorrVals(corrSets,opt1,opt2)

% check if the given corrSets are cells or edges:
if isfield(corrSets,'resF')
    isCells=1;
else
    isCells=0;
end

if isCells
    for idx=1:length(corrSets)
        flags=horzcat(corrSets(idx).edge.flag);
        checkVec=(sum(flags,2)==length(corrSets(idx).edge));
        for edgeId=1:length(corrSets(idx).edge)
            if strcmp(opt1,'usefm')
                corrSets(idx).edge(edgeId).fcorr=corrSets(idx).edge(edgeId).fm;
            elseif strcmp(opt1,'usefn')
                corrSets(idx).edge(edgeId).fcorr=NaN+zeros(length(checkVec),2);
                corrSets(idx).edge(edgeId).fcorr(checkVec,:)=corrSets(idx).edge(edgeId).fm(checkVec,:);
            elseif strcmp(opt1,'usefc')
                corrSets(idx).edge(edgeId).fcorr=corrSets(idx).edge(edgeId).fc;
            else
                error('The given option is not supported');
            end
        end
    end
else
    for idx=1:length(corrSets)
        flags=horzcat(corrSets(idx).flag);
        checkVec=(flags==1);
        
        if strcmp(opt1,'usefm')
            corrSets(idx).fcorr=corrSets(idx).fmMag;
        elseif strcmp(opt1,'usefn')
            corrSets(idx).fcorr=NaN+zeros(length(checkVec),1);
            corrSets(idx).fcorr(checkVec,1)=corrSets(idx).fmMag(checkVec);
        elseif strcmp(opt1,'usefc')
            corrSets(idx).fcorr=corrSets(idx).fcMag;
        else
            error('The given option is not supported');
        end
        
        if strcmp(opt2,'useItot')
            corrSets(idx).Icorr=corrSets(idx).Itot;
        elseif strcmp(opt2,'useIavg')
            corrSets(idx).Icorr=corrSets(idx).Iavg;
        else
            error('The given option is not supported');
        end
    end
end