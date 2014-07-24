function varargout=collectEdgeValues(groupedClusters,goodSet,varargin)
% [deg_vals,lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'deg','lgth','f1','f2','fc1','fc2')
% Runs through the groupedClusters and collects all the data from fields
% specified in the input arguments. Potential fields are:
% 'deg'       : The degree of connectivity of a cell.
% 'lgth'      : The interface length.
% 'dlpix'     : Returns the distance between stress-measurement points
%               along the whole interface (not only internal).
% 'f1'        : The force f1 determined by the network analysis
% 'f2'        : The force f2 determined by the network analysis
% 'fc1'       : The force f1 determined by the network analysis
% 'fc2'       : The force f2 determined by the network analysis
% 'nVec'      : is actually nVec_internal, the normal on the internal part of
%               the interface (and not the field n_Vec which contains segments
%               outside the cell).
% 'Itot'      : The total Ecad intensity along the interface (raw image average
%               with average filter of radius r1=10pix)
% 'Iavg'      : The total Ecad intensity along the interface (raw image average
%               with average filter of radius r1=10pix)
% 'SIcorr'    : The stress and intensity profile along the INNER interface.
%               Spatial but no temporal information.
% 'SIcorrRand': The stress and intensity profile along the INNER interface.
%               Needs two more inputs:
%                  1) the length of the interface over which the values 
%                     should be randomized. 
%                  2) if the values should be randomized or if just the 
%                     regrouping should be performed. This is helpful for
%                     comparing the randomized correlation coefficients
%                     with original pairings on exactly the same data.
%                     (note that subdividing the interface into smaller
%                     pieces requires either abondoning are doublecounting
%                     data.)
%               Spatial but no temporal information.
% 'intfRand'  : randomizes intensity-stress pairs over one interface in a
%               specific frame. 
% 'frameRand' : randomizes intensity-stress pairs over all interfaces
%               within one cluster at one frame. 
% 'fullRand'  : randomizes all intensity-stress pairs irrespective of
%               interface, cluster or frame, i.e. completely random
%               pairings.
% 'corr'      : structure for cross-correlating intensity and
%               cell-cell-forces over time.
% 'r2'        : if option 'r2' is given we take r2=5pix values for averaging
%               the Ecad images. This applies to 'Itot', 'Iavg', 'SIcorr'
%               and 'SIcorrRand
% 'IvarMin'   : if option 'IvarMin' is given, we consider only edges whose
%               intensity variations are above a certain cut-off value.
%               this cut-off value needs to be specified as second
%               arguement. Acts on SIcorr, SIcorrRand, intfRand, frameRand,
%               fullRand. Returns the values for the intensity variances
%               that have been used.
% 'SIcorrCoef': returns the correlation coefficients for stress and
%               intensity at each interface, the p-value, upper and lower 
%               bounds for the confidence intervall, as well as the the max
%               difference between upper/lower bound and the correlation
%               coefficient. Values are calculated for each frame.
% 'SIcorrCoefIntfRand': Randomizes all pairs of an interface in a specific
%               frame and returns the correlation coefficient as
%               SIcorrCoef.
% 'SIcorrCoefRand': works as 'SIcorrRand' but returns the correlation 
%               coefficients for stress and intensity at each interface,
%               see 'SIcorrCoef'. Values are calculated at each frame.


degPos=find(strcmp('deg',varargin));
if ~isempty(degPos)
    degCheck = 1;
    deg_vals = [];
else
    degCheck = 0;
end
    

lgthPos=find(strcmp('lgth',varargin));
if ~isempty(lgthPos)
    lgthCheck = 1;
    % it is the next entry which contains the numeric value:
    lgth_vals   = [];
else
    lgthCheck = 0;
end

dlpixPos=find(strcmp('dlpix',varargin));
if ~isempty(dlpixPos)
    dlpixCheck = 1;
    % it is the next entry which contains the numeric value:
    dlpix_vals   = [];
else
    dlpixCheck = 0;
end


f1Pos=find(strcmp('f1',varargin));
if ~isempty(f1Pos)
    f1Check = 1;
    % it is the next entry which contains the numeric value:
    f1_vals   = [];
else
    f1Check = 0;
end

f2Pos=find(strcmp('f2',varargin));
if ~isempty(f2Pos)
    f2Check = 1;
    % it is the next entry which contains the numeric value:
    f2_vals   = [];
else
    f2Check = 0;
end

fc1Pos=find(strcmp('fc1',varargin));
if ~isempty(fc1Pos)
    fc1Check = 1;
    % it is the next entry which contains the numeric value:
    fc1_vals   = [];
else
    fc1Check = 0;
end

fc2Pos=find(strcmp('fc2',varargin));
if ~isempty(fc2Pos)
    fc2Check = 1;
    % it is the next entry which contains the numeric value:
    fc2_vals   = [];
else
    fc2Check = 0;
end

nVecPos=find(strcmp('nVec',varargin));
if ~isempty(nVecPos)
    nVecCheck = 1;
    % it is the next entry which contains the numeric value:
    nVec_vals   = [];
else
    nVecCheck = 0;
end

ItotPos=find(strcmp('Itot',varargin));
if ~isempty(ItotPos)
    ItotCheck = 1;
    % it is the next entry which contains the numeric value:
    Itot_vals   = [];
else
    ItotCheck = 0;
end

IavgPos=find(strcmp('Iavg',varargin));
if ~isempty(IavgPos)
    IavgCheck = 1;
    % it is the next entry which contains the numeric value:
    Iavg_vals   = [];
else
    IavgCheck = 0;
end

SIcorrPos=find(strcmp('SIcorr',varargin));
if ~isempty(SIcorrPos)
    SIcorrCheck = 1;
    % it is the next entry which contains the numeric value:
    SIcorr_vals   = [];
else
    SIcorrCheck = 0;
end

SIcorrCoefPos=find(strcmp('SIcorrCoef',varargin));
if ~isempty(SIcorrCoefPos)
    SIcorrCoefCheck = 1;
    SIcorrCoef_vals = [];
    SIcorr_vals     = [];
else
    SIcorrCoefCheck = 0;
end

SIcorrCoefIntfRandPos=find(strcmp('SIcorrCoefIntfRand',varargin));
if ~isempty(SIcorrCoefIntfRandPos)
    SIcorrCoefIntfRandCheck = 1;
    SIcorrCoef_vals = [];
    SIcorr_vals     = [];
else
    SIcorrCoefIntfRandCheck = 0;
end

if SIcorrCoefCheck && SIcorrCoefIntfRandCheck
    error('Options SIcorrCoef and SIcorrCoefIntfRand cannot be used together')
end

SIcorrCoefRandPos=find(strcmp('SIcorrCoefRand',varargin));
if ~isempty(SIcorrCoefRandPos)
    SIcorrCoefRandCheck = 1;
    % it is the next entry which contains the numeric value:
    SIcorrRand_lgth        = varargin{SIcorrCoefRandPos+1};
    SIcorrRand_onlyRegroup = varargin{SIcorrCoefRandPos+2};
    SIcorrCoefRand_vals = [];
    SIcorrRand_vals = [];
    SIcorr_vals     = [];
else
    SIcorrCoefRandCheck = 0;
end

SIcorrRandPos=find(strcmp('SIcorrRand',varargin));
if ~isempty(SIcorrRandPos)
    SIcorrRandCheck = 1;
    % it is the next entry which contains the numeric value:
    SIcorrRand_lgth        = varargin{SIcorrRandPos+1};
    SIcorrRand_onlyRegroup = varargin{SIcorrRandPos+2};
    SIcorrRand_vals = [];
    SIcorr_vals     = [];
else
    SIcorrRandCheck = 0;
end

intfRandPos=find(strcmp('intfRand',varargin));
if ~isempty(intfRandPos)
    intfRandCheck = 1;
    % the correlation pairs have to be drawn from the data structure
    SIcorr_vals     = [];
    SIcorr_vals_intfRand = [];
else
    intfRandCheck = 0;
end


frameRandPos=find(strcmp('frameRand',varargin));
if ~isempty(frameRandPos)
    frameRandCheck = 1;
    % the correlation pairs have to be drawn from the data structure
    SIcorr_vals     = [];
    SIcorr_vals_ext = [];
    maxClusterId=-Inf;
else
    frameRandCheck = 0;
end

fullRandPos=find(strcmp('fullRand',varargin));
if ~isempty(fullRandPos)
    fullRandCheck = 1;
    % the correlation pairs have to be drawn from the data structure
    SIcorr_vals   = [];
else
    fullRandCheck = 0;
end

r2Pos=find(strcmp('r2',varargin));
if ~isempty(r2Pos)
    r2Check = 1;
else
    r2Check = 0;
end

IvarMinPos=find(strcmp('IvarMin',varargin));
if ~isempty(IvarMinPos)
    IvarMinCheck = 1;
    IvarMin      = varargin{IvarMinPos+1};
    Ivar_vals    = [];
else
    IvarMinCheck = 0;
end

corrPos=find(strcmp('corr',varargin));
if ~isempty(corrPos)
    corrCheck = 1;
    % initialize:
    maxIdx =length(goodSet);
    corr_out(maxIdx).fcMag= [];
    corr_out(maxIdx).fmMag= [];
    corr_out(maxIdx).Itot = [];
    corr_out(maxIdx).Iavg = [];
    corr_out(maxIdx).t    = [];
else
    corrCheck = 0;
end

for idx=1:length(goodSet)
    clusterId=goodSet(idx).clusterId;
    edgeId   =goodSet(idx).edgeId;
    toDoList =goodSet(idx).frames;
    
    % initialize all values:
    corr_out(idx).clusterId = clusterId;
    corr_out(idx).edgeId    = edgeId;
    corr_out(idx).t         = NaN+zeros(toDoList(end),1);
    corr_out(idx).fcMag     = NaN+zeros(toDoList(end),1);
    corr_out(idx).fmMag     = NaN+zeros(toDoList(end),1);
    corr_out(idx).Itot      = NaN+zeros(toDoList(end),1);
    corr_out(idx).Iavg      = NaN+zeros(toDoList(end),1);
    corr_out(idx).flag      = ones(toDoList(end),1);
    corr_out(idx).frames    = 1:toDoList(end);
    

    
    for frame=toDoList
        % Now go through all checks:
        if  degCheck
            nodes=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.nodes;
            deg_vals=vertcat(deg_vals,[groupedClusters.cluster{clusterId}.trackedNet{frame}.node{nodes(1)}.deg groupedClusters.cluster{clusterId}.trackedNet{frame}.node{nodes(2)}.deg]);
        end
        
        if  lgthCheck
            lgth_vals=vertcat(lgth_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.intf_internal_L);
        end
        
        if  dlpixCheck
            % full_length = groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.intf;
            curr_cntrs=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.cntrs;
            curr_dlpix=sqrt(sum((curr_cntrs(2:end,:)-curr_cntrs(1:end-1,:)).^2,2));
            dlpix_vals=vertcat(dlpix_vals ,curr_dlpix);
        end
        
        if  f1Check
            f1_vals=vertcat(f1_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f1);
        end
        
        if  f2Check
            f2_vals=vertcat(f2_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f2);
        end
        
        if  fc1Check
            fc1_vals=vertcat(fc1_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.fc1);
        end
        
        if  fc2Check
            fc2_vals=vertcat(fc2_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.fc2);
        end
        
        if  nVecCheck
            nVec_vals=vertcat(nVec_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.nVec_internal);
        end
        
        if  ItotCheck
            if r2Check
                Itot_vals=vertcat(Itot_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.tot2);
            else % this is the default:
                Itot_vals=vertcat(Itot_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.tot);
            end
        end
        
        if  IavgCheck
            if r2Check
                Iavg_vals=vertcat(Iavg_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.avg2);
            else
                Iavg_vals=vertcat(Iavg_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.avg);
            end
        end
        
        % this is the correlation of the proper pairing
        if  SIcorrCheck || SIcorrRandCheck || fullRandCheck || frameRandCheck || intfRandCheck || SIcorrCoefCheck || SIcorrCoefIntfRandCheck || SIcorrCoefRandCheck
            SPosIn =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.cntrs;
            SValsIn=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.s_vec;
            IPosIn =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.intf_internal;
            if r2Check
                if isfield(groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int,'val2') 
                    IValsIn=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.val2;
                else
                    display('had to omit data as val2 field was missing:')
                    display(['setId: ' ,num2str(idx),';   frame: ',num2str(frame)])
                    IValsIn=[]; SValsIn=[]; SPosIn=[]; IPosIn=[];
                end
            else
                IValsIn=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.val;
            end
            [IVals,~,SMagVals,~]=pairSIprofiles(SPosIn,SValsIn,IPosIn,IValsIn);
            
            if IvarMinCheck
                %calculate the variance of the intensity along the
                %interface:
                curr_Ivar=var(IVals);
                if curr_Ivar<IvarMin
                    % dump the values:
                    SMagVals=[]; IVals=[];
                    display('data omitted:')
                    display(['unsufficient Ivar, cluster: ',num2str(clusterId),'; edge: ',num2str(edgeId),'; frame: ',num2str(frame)])
                else
                    % data set is used as is, store the used variance:
                    Ivar_vals=vertcat(Ivar_vals ,curr_Ivar);
                end
            end
            
            % load the values if there are still any (IvarMinCheck):
            SIcorr_vals=vertcat(SIcorr_vals ,[SMagVals IVals]);
            
            % only if intfRandCheck is given:  
            if intfRandCheck
                numVals=length(IVals);
                messVec = randperm(numVals);
                % of course we change the order only of one vec!
                IVals_rand    = IVals(messVec);
                % whereas SMagVals will always keep the old order:
                SIcorr_vals_intfRand=vertcat(SIcorr_vals_intfRand ,[SMagVals IVals_rand]);
            end
            
            % only if frameRandCheck is given:    
            if frameRandCheck
                numVals=length(IVals);
                initvec=zeros(numVals,1);
                clusterIDvec=clusterId+initvec;
                frameIDvec  =frame+initvec;
                SIcorr_vals_ext=vertcat(SIcorr_vals_ext ,[SMagVals IVals clusterIDvec frameIDvec]);
                % will be permutated at the very end, see below!
            end            
        end
        
        if (SIcorrCoefCheck || SIcorrCoefIntfRandCheck) && ~isempty(SMagVals)
            %calculate the numerical correlation coefficients control:
            if SIcorrCoefIntfRandCheck
                numSMagVals=length(SMagVals);
                messVec = randperm(numSMagVals);
                % of course we change the order only of one vec!
                SMagVals_rand    = SMagVals(messVec);
            else
                SMagVals_rand    = SMagVals;
            end
            [r_SI p_SI rLO_SI rUP_SI] = corrcoef(SMagVals_rand,IVals);
            r_SI_STD=max(max([r_SI-rLO_SI rUP_SI-r_SI]));
            
            % also append the clusterId, edgeId to later on determine from
            % how many different edegs and clusters the data has been drawn.
            SIcorrCoef_vals=vertcat(SIcorrCoef_vals ,[r_SI(1,2) p_SI(1,2) rLO_SI(1,2) rUP_SI(1,2) r_SI_STD clusterId edgeId]);
        end
        
        % this is the correlation of the messed-up pairings
        if  SIcorrRandCheck || SIcorrCoefRandCheck
            % the proper pairings have been established above. Now we mess
            % them up as wanted:
            cum_pairs=[];
            numIVals=length(IVals);
            for istPt=1:SIcorrRand_lgth:numIVals
                % istPt is the start point of the current interval, the
                % intervall is SIcorrRand_lgth entries long. 
                % clear the values to make sure that old values are not used multiple times.
                curr_IVals=[]; 
                curr_SMagVals=[];
                curr_IVals_rand=[];
                curr_SMagVals_rand=[];
                
                % Pull the current values:
                endPt=istPt+SIcorrRand_lgth-1;
                if endPt<=numIVals
                    curr_IVals    = IVals(istPt:endPt);
                    curr_SMagVals = SMagVals(istPt:endPt);
                elseif numIVals>=SIcorrRand_lgth
                    % the last interval was shorter than SIcorrRand_lgth,
                    % so we flip the interval:
                    curr_IVals    = IVals(end-SIcorrRand_lgth+1:end);
                    curr_SMagVals = SMagVals(end-SIcorrRand_lgth+1:end);
                else
                    display('interface had to be dismissed as it is shorter than SIcorrRand_lgth')
                    display(['setId: ' ,num2str(idx),';   frame: ',num2str(frame)])
                end
                % now messup the order of the current value pairs (if we found enough):
                if ~isempty(curr_IVals) %==SIcorrRand_lgth
                    messVec = randperm(SIcorrRand_lgth);
                    % of course we change the order only of one vec!
                    if SIcorrRand_onlyRegroup
                        % this must be used to get unrandomized parings:
                        curr_IVals_rand    = curr_IVals;
                    else
                        % here, IVals will be mixed up:
                        curr_IVals_rand    = curr_IVals(messVec);
                    end
                    % whereas SMagVals will always keep the old order:
                    curr_SMagVals_rand = curr_SMagVals;
                    % taken out:
                    % SIcorrRand_vals=vertcat(SIcorrRand_vals ,[curr_SMagVals_rand curr_IVals_rand]);
                    % newline:
                    cum_pairs=vertcat(cum_pairs ,[curr_SMagVals_rand curr_IVals_rand]);
                end
            end
            % newline:
            SIcorrRand_vals=vertcat(SIcorrRand_vals ,cum_pairs);
        end

        % uses the code from above!!!
        if SIcorrCoefRandCheck && ~isempty(cum_pairs)
            % calculate the numerical correlation coefficients.
            % Randomization (or not) is managed above in the main
            % section of SIcorrRandCheck):
            [r_SI p_SI rLO_SI rUP_SI] = corrcoef(cum_pairs(:,1),cum_pairs(:,2));
            r_SI_STD=max(max([r_SI-rLO_SI rUP_SI-r_SI]));
            
            % also append the clusterId, edgeId to later on determine from
            % how many different edegs and clusters the data has been drawn.
            SIcorrCoefRand_vals=vertcat(SIcorrCoefRand_vals ,[r_SI(1,2) p_SI(1,2) rLO_SI(1,2) rUP_SI(1,2) r_SI_STD clusterId edgeId]);
        end
        
        if  corrCheck
            % collect the intensity values:
            if r2Check
                corr_out(idx).Itot(frame,1) =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.tot2;
                corr_out(idx).Iavg(frame,1) =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.avg2;
            else % this is the default:
                corr_out(idx).Itot(frame,1) =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.tot;
                corr_out(idx).Iavg(frame,1) =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.avg;
            end
            
            % Determine fmMag and fcMag. The direction is not important!
            f1 =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f1;
            f2 =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f2;
            fc =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.fc;            
            fn = 0.5*(f1-f2);
            
            corr_out(idx).fcMag(frame,1)=norm(fc);            
            % create the mixed force vector which takes the network
            % force whenever possible.
            if ~isnan(sum(fn))
                corr_out(idx).fmMag(frame,1) = norm(fn);
                % to keep track which value we have chosen:
                corr_out(idx).flag(frame,1)=1;
            else
                corr_out(idx).fmMag(frame,1) = norm(fc);
                % to keep track which value we have chosen:
                if ~isnan(sum(fc))
                    corr_out(idx).flag(frame,1)=0;
                else
                    % then, the edge was empty and it is OK to use the
                    % network result
                    corr_out(idx).flag(frame,1)=1;
                end
            end
            % The value will be overwritten many times but in this way we
            % don't have to check which is the first non-empty entry in the
            % trackedNet structure.
            corr_out(idx).dt_mean = groupedClusters.cluster{clusterId}.trackedNet{frame}.par.dt_mean;
            corr_out(idx).dt_std  = groupedClusters.cluster{clusterId}.trackedNet{frame}.par.dt_std;
            
            % read out also the absolute time point:
            corr_out(idx).t(frame,1)=groupedClusters.cluster{clusterId}.trackedNet{frame}.par.t;
        end
    end
end

if fullRandCheck
    numPairs=length(SIcorr_vals(:,1));
    messVec = randperm(numPairs);
    % messup one column of the pairings:
    SIcorr_vals(:,1)=SIcorr_vals(messVec,1);
end

if frameRandCheck
    maxClusterId=max(SIcorr_vals_ext(:,3));
    maxFrameId  =max(SIcorr_vals_ext(:,4));
    for iCluster=1:maxClusterId
        for iFrame=1:maxFrameId            
            cv1=(SIcorr_vals_ext(:,3)==iCluster) & (SIcorr_vals_ext(:,4)==iFrame);
            if ~isempty(SIcorr_vals_ext(cv1,:))
                % display(['cluster id: ',num2str(iCluster),';    frame id: ',num2str(iFrame)])
                % these are the pair-IDs that belong to the same cluster
                % and frame:
                currPairIds=find(cv1);
                numPairs=length(currPairIds);
                messVec = randperm(numPairs);
                %messup the order of the pairIds
                currPairIds_messed=currPairIds(messVec);
                % messup the current pairs by associating randomized stress 
                % with original intensities values (from this subgroup
                % only):
                SIcorr_vals_ext(currPairIds,1:2)=[SIcorr_vals_ext(currPairIds_messed,1) SIcorr_vals_ext(currPairIds,2)];
            end
        end
    end
end

% number of optional input = number of output arguments:
optargin = size(varargin,2);


if degCheck
    varargout(degPos)  = {deg_vals};
end

if lgthCheck
    varargout(lgthPos)  = {lgth_vals};
end

if dlpixCheck
    varargout(dlpixPos)  = {dlpix_vals};
end

if f1Check
    varargout(f1Pos) = {f1_vals};
end

if f2Check
    varargout(f2Pos) = {f2_vals};
end

if fc1Check
    varargout(fc1Pos) = {fc1_vals};
end

if fc2Check
    varargout(fc2Pos) = {fc2_vals};
end

if nVecCheck
    varargout(nVecPos) = {nVec_vals};
end

if ItotCheck
    varargout(ItotPos) = {Itot_vals};
end

if IavgCheck
    varargout(IavgPos) = {Iavg_vals};
end

if SIcorrCheck
    varargout(SIcorrPos) = {SIcorr_vals};
end

if SIcorrCoefCheck || SIcorrCoefIntfRandCheck
    varargout(SIcorrCoefPos)   = {SIcorrCoef_vals};
end

if SIcorrCoefIntfRandCheck
    varargout(SIcorrCoefIntfRandPos)   = {SIcorrCoef_vals};
end

if SIcorrRandCheck
    varargout(SIcorrRandPos)   = {SIcorrRand_vals};
    varargout(SIcorrRandPos+1) = {SIcorrRand_lgth};
    varargout(SIcorrRandPos+2) = {SIcorrRand_onlyRegroup};
end

if SIcorrCoefRandCheck
    varargout(SIcorrCoefRandPos)   = {SIcorrCoefRand_vals};
    varargout(SIcorrCoefRandPos+1) = {SIcorrRand_lgth};
    varargout(SIcorrCoefRandPos+2) = {SIcorrRand_onlyRegroup};
end

if fullRandCheck
    varargout(fullRandPos) = {SIcorr_vals};
end

if frameRandCheck
    varargout(frameRandPos) = {SIcorr_vals_ext(:,1:2)};
end

if intfRandCheck
    varargout(intfRandPos) = {SIcorr_vals_intfRand};
end

if r2Check
    r2_vals='taken r2';
    varargout(r2Pos) = {r2_vals};
end

if IvarMinCheck
    varargout(IvarMinPos)   = {Ivar_vals};
    varargout(IvarMinPos+1) = {Ivar_vals};
end

if  corrCheck
    varargout(corrPos) = {corr_out};
end

if length(varargout)~=optargin
    error('Something went wrong');
end

for k=2:optargin
    if size(varargout{1},1)~=size(varargout{k},1)
        display('Something might have went wrong, vectors are of different length!');
    end
end