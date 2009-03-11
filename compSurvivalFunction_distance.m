function [survFunc, pvec, decaytimes] = compSurvivalFunction_distance(data1, censor, distancevec);
% compSurvivalFunction_general compares the survival functions for
% different conditions
%
% SYNOPSIS [survFunc,pvec] = compSurvivalFunction_general(data1, field1,
% data2, field2)
%
% INPUT     data        = experiment structure with data set
%           censor      = censor variable
%                       if censor==0 use all trajectories, both the complete
%                       and the cut off ones
%                       if censor==1 (default) use only the complete
%                       trajectories
%           distancevec = distance vector, e.g. [-10:10:50]
%                       specifies the distance bins (in pixel) from the
%                       segmented edge that are considered separately
%                       NOTE: negative distances refer to the areas inside
%                       the pattern, positive distance to the distance
%                       outside the edge
%                   
%                   
% OUTPUT    survFunc    = survival function matrix with survival functions
%                       for the different distance bins                    
%           pvec        = matrix with p-values determined from a KS-test of
%                       the survival functions of each distance bin with
%                       every other distance bin
%                     
%           decaytimes  = decay times, i.e. number of frames for which the
%                       respective survival function reaches the percentage
%                       levels specified by percvec; this has the format
%                   
% NOTE: This function works best if the distances from the segmented have
% been determined previosuly and stored in the field data.segmentEUdist; if
% this field does not yet exist, run the function 
% data = fillStructSegmentStatus(data1);
% to calculate them.
%
%
% last modified DATE: 28-Aug-2008 (Dinah)
% last modified DATE: 10-Sep-2008 (Dinah)
% last modified DATE: 02-Oct-2008 (Dinah)
% last modified DATE: 04-Nov-2008 (Dinah)


% averaging can only be performed up until the minimum common length of all 
% movies, so we determine the shortest movie length in the structure
% averaging can only be performed up until the minimum common length of all 
% movies, so we determine the shortest movie length in the structure
n1 = length(data1);

for i=1:n1
    mlvec1(i) = data1(i).movieLength;
end
minlen = min(mlvec1);


censoring = 0;
if nargin>1
    if ~isempty(censor)
        censoring = censor;
    end
end

detEU = 1;
if isfield(data1,'segmentEUdist')
    if ~isempty(data1(1).segmentEUdist)
        detEU = 1;
    end
end

if detEU==0
    disp('Euclidian distance vector do not exist yet, so the function is calculating them now...');
    data1 = fillStructSegmentStatus(data1);
end


od = cd;

for i=1:length(data1)
    
    fprintf('movie #%02d',i);
        
    
    %%=====================================================================
    % read lifetime and lifetime status matrices from LifetimeInfo
    %%=====================================================================
       
    % current path
    path = data1(i).source;

    % number of frames for this exp
    lenf = data1(i).movieLength;

    % load lifetime inof data file
    lftpath = [path,'/LifetimeInfo'];
    cd(lftpath);
    lftname = 'lftInfo.mat';
    loadfile = load(lftname);
    cd(od);

    lftInfo = loadfile.lftInfo;
    
    
    %%=====================================================================
    % read (or determine, if necessary) pattern segmentation status
    %%=====================================================================
    
    % check if segmentation status already exists
    if isfield(data1,'segmentStatus')
        
        segmentStatusVector = data1(i).segmentStatus;
    % else determine the segmentation status here and fill in the value    
    else
              
        % load segmentation image from specified location
        SegmFileName = data1(i).segmentDataFileName;
        SegmFilePath = data1(i).segmentDataFilePath;

        cd(SegmFilePath);
        SegmentMask = imread(SegmFileName);    

        cd(od);

        % calculate segmentation status (1=Inside, 0=oustide segmented region)
        [segmentStatusVector] = calcIORegionLfthistSimple(lftInfo, SegmentMask);
        
        data1(i).segmentStatus = segmentStatusVector';
        
    end
    
    %%=====================================================================
    % read pattern segmentation Euclidian distance vector
    %%=====================================================================
    
    segmentDistanceVector = data1(i).segmentEUdist;
    
    
    %%=====================================================================
    % read out desired lifetimes
    %%=====================================================================
    
    lftMat = lftInfo.Mat_lifetime;
    statMat =  lftInfo.Mat_status;
    [sx,sy] = size(lftMat);
    lftVec = nan*zeros(sx,1);
    
    % IF censoring==1: a trajectory is counted for the lifetime analysis if
    % the status of the trajectory is ==1 AND the value of any gaps is ==4
    % IF censoring==0: count all trajectories of all status values (1,2,3)
    % while the gap values are ==4
    for k=1:sx
        % current status vector
        cstat = nonzeros(statMat(k,:));
        % current lifetime vector
        clft = lftMat(k,:);
        
        % counting status (count or don't count this entry)
        countStat = ( (min(cstat)==1) & (max(cstat)<5) );
        if censoring==0
            countStat = (max(cstat)<5);
        end
        
        if countStat==1
            lftVec(k) = max(clft);
        end
            
    end % of for k-loop    
    
    % lifetime vectors containing all counted lifetime lengths
    SegmentedLFTInfoMat = nan*zeros(sx,3);
    SegmentedLFTInfoMat(:,1) = lftVec;
    SegmentedLFTInfoMat(:,2) = segmentStatusVector;
    SegmentedLFTInfoMat(:,3) = max(segmentDistanceVector,[],2);

    % lifetime histograms
    data1(i).segmentedLifetimeInfo = SegmentedLFTInfoMat;
    
    fprintf('\b\b\b\b\b\b\b\b\b');  
    
    if i==1
        segmentedLFTmatrix = SegmentedLFTInfoMat;
    else
        segmentedLFTmatrix = [segmentedLFTmatrix; SegmentedLFTInfoMat];
    end

end


% find inside pattern positions
fipos = find(segmentedLFTmatrix(:,2) == 1);
% set distances to negative
segmentedLFTmatrix(fipos,3) = - segmentedLFTmatrix(fipos,3);
% sort matrix
segmentedLFTsort = sortrows(segmentedLFTmatrix,3);
segmentedLFTcurr = segmentedLFTsort;

for n=1:length(distancevec)-1
    cdist = segmentedLFTcurr(:,3);
    cpos = find( (cdist>= distancevec(n)) & (cdist<distancevec(n+1)) );
    clft = segmentedLFTcurr(cpos,1);
    chist = hist(clft,[1:minlen]);
    histmat(n,:) = chist;
    segmentedLFTcurr(1:max(cpos),:)=[];
end

survFuncMat = repmat(sum(histmat,2),1,minlen) - cumsum(histmat,2);


%% p-vales
% convert to distributions
for i=1:length(distancevec)-1
    df1 = convSurvivalFunction2dist(survFuncMat(i,:));
    for k=i:length(distancevec)-1
        df2 = convSurvivalFunction2dist(survFuncMat(k,:));
        [H,pval_ik] = kstest2(df1,df2);
        pvalMat(i,k) = pval_ik;
    end
end

histmat1 = survFuncMat;

%% level output

levels = [0.9:-0.1:0.1];
if nargin>6
    if ~isempty(percvec)
        levels = percvec;
    end
end

for k=1:length(levels)
        
    % current percentage increment
    cinc = levels(k);
    % find decay time for this level 
    for a=1:length(distancevec)-1
        fpos = find( (histmat1(a,:)/histmat1(a,1))<cinc );
        if ~isempty(fpos)
            incrVec1(a,k) = min( fpos );
        else
            incrVec1(a,k) = nan;
        end                  
    end
    
end

%% ==========================================================================
%% display results


figure; hold on;
for n=1:length(distancevec)-1
    histmat1norm(n,:) = histmat1(n,:)/histmat1(n,1);
end
plot(histmat1norm');
xlabel('lifetime');
ylabel('survival fct');


figure; hold on;
Mat1 = incrVec1;
Mat2 = repmat([distancevec(1:length(distancevec)-1)]',1,length(levels));
plot(Mat2,Mat1);
xlabel('distance in pixel from pattern edge');
ylabel('lifetime of level');


survFunc = histmat1norm;
decaytimes = incrVec1;
pvec = pvalMat;


end % of function




%=========================================================================
%
%                            SUBFUNCTION



function [df] = convSurvivalFunction2dist(sf)
% convSurvivalFunction2dist converts a survival function back into a
% distribution function
%
% SYNOPSIS [df] = convSurvivalFunction2dist(sf)
%
% INPUT     sf    = survival function, presumed at equal spacing 
%
% OUTPUT    df    = distribution function
%
% NOTE: Advantage of this function is that survival functions can be summed
% and cropped until appropriate time points, that that the reconversion
% gets rid of the problem of unequal movie lengths


% convert to histogram
hf = abs(diff(sf));

df = zeros(sum(hf),1);
ct = 0;

for i=1:length(hf)
    numentry = round(hf(i));
    %cv = i + zeros(numentry,1);
    df(ct+1:ct+numentry) = i; 
    ct = ct+numentry;
end % of for

end % of subfunction




    
