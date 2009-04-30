function [results] = stageFitLifetimes_E2E(data)
% stageFitLifetimes_E2E fits the lifetimes contained in the structure in 
% several stages, using the slow and fast data, like the function
% stageFitLifetimes - in this function, individual movies from the data
% structure are removed one by one (providing there's enough movies), in
% order to allow a jackknife error estimation in a subsequent step
%
% SYNOPSIS [results] = stageFitLifetimes_E2E(data)
%
% INPUT     data:   structure containing all raw data
%           
%
% OUTPUT:   results: structure containing the fields
%                       .hist_fast
%                       .hist_slow
%                       .numcells_slow 
%                       .numcells_fast
%           lftHist_fast =  2xn vector containing 
%                           row 1: time vector
%                           row 2: normalized lifetimes 
%                           (for FAST acquisition)
%           lftHist_slow =  2xn vector containing 
%                           row 1: time vector
%                           row 2: normalized lifetimes 
%                           (for SLOW acquisition)
%           NOTE: individual structure fields contain the results
%           corresponding to the removal of subsequent movies from the data
%
% last modified DATE: 31-Jul-2007 (Dinah)


%% ========================================================================
% 
%       FIRST STEP: determine which entries in the data structure represent
%       fast or slow movies, so that we know which ones of the movies can
%       be deleted
%
% =========================================================================


% % initialize variables
% ct_fast = 1;
% ct_slow = 1;
% 
% % In this implementation, since we want to average over the data without
% % any breaks due to averaging over different numbers of histograms, the
% % histograms are all cut to the same length. Thus, we need to know what's 
% % the maximum common time span for all the fast/slow movies - as we only
% % input movies with the same framerate (2s slow, 0.4s fast), the knowledge
% % of the last frame suffices
% 
% % loop over all entries in the data structure
% for i=1:length(data)
%     
%     % read current detection frequency
%     framerate_i = data(i).framerate;
%         
%     %  if movie is fast and data exists
%     if ( (framerate_i==0.4) & ~isempty(data(i).lftHist) )
%         % enter i as a position of a fast movie
%         positions_fast(ct_fast) = i;
%         ct_fast = ct_fast+1;
%         
%     %  if movie is slow and data exists    
%     elseif ( (framerate_i==2) & ~isempty(data(i).lftHist) )
%         % enter i as a position of a slow movie
%         positions_slow(ct_slow) = i;
%         ct_slow = ct_slow+1;  
%          
%     end % of if
% end
% 
% 
% ct_del = 1;
% 
% for i=1:length(data)
%     
%     % read current detection frequency
%     framerate_i = data(i).framerate;
%     
%     
%     % if more than one movie for this condition (fast or slow) exists,
%     % then delete this movie and perform data analysis without it
%     
%     deletevar = 0;
%     % test if more than one movie of specified condition exists
%     if ( (framerate_i==2) & (ct_slow>1) )
%         deletevar = 1;
%     elseif ( (framerate_i==0.4) & (ct_fast>1) )
%         deletevar = 1;
%     end
% 
%     % if so, delete movie i and analyse the rest    
%     if deletevar == 1
%         
%         data_current = data;
%         data_current(i) = [];
%         
%         [results_current] = stageFitLifetimesPlat(data_current);
%         results(ct_del) = results_current;
%         ct_del = ct_del+1;
%     end
%     
% end % of for



%% ========================================================================
%  persistence cutoff: all pits that live longer than pc are considered
%  persistent (meas in seconds)

%pc = persThresh;


%% ========================================================================
% 
%       FIRST STAGE: average all data to identify initial detection
%       artifact 
%
% =========================================================================

xresample = 0;

cutRes = 0;
if nargin>1
    if xresample==1
        cutRes = 1;
    end
end
        

% NOTE: since we're interested in the initial decay of the detection
% artifact, it should be sufficient to fit only the 100 or so frames; the
% cutoff is defined as cutoff_art

cutoff_art = 100;

% initialize results matrix
histMatrix_all = zeros(length(data),cutoff_art);

for i=1:length(data)
    
    % read lifetime Histogram
    currHist = data(i).lftHist;
    
    if ~isempty(currHist)            
        % normalize over area betwen 0 and cutoff_art frames
        currHistNorm = currHist(1:cutoff_art)/sum(currHist(1:cutoff_art));
        % enter into matrix of all results
        histMatrix_all(i,1:length(currHistNorm)) = currHistNorm;
    end
end

% average results...
histVector_all = nanmean(histMatrix_all,1);
%...and renormalize
histVectorRenorm_all = histVector_all/nansum(histVector_all);

tvec = [1:cutoff_art];



% plot results
f1 = figure;
plot(tvec,histVectorRenorm_all,'b.-'); hold on;
xlabel('frames'); ylabel('average frequency');
axis([0 cutoff_art 0 1.1*max(histVectorRenorm_all(:))]);



% fit averaged results with multiple exponentials (weibulls with k==1)
guessvector = [0    0.3 0.3 1   0.3 5   1   0.3 20  1];
fixvector =   [0    0   0   1   0   0   1   0   0   1]; 
[estAll, resAll] = fitcurveMultiWeibullODF_lsq(tvec, histVectorRenorm_all, guessvector, fixvector);
% take abs of all except offset
estAll(2:length(estAll)) = abs(estAll(2:length(estAll)));


% plot first component, which is the detection artifact
currEstArt = abs(estAll(2))*(1/estAll(3))*exp(-tvec/estAll(3));
hold on; plot(tvec,currEstArt,'g-');
title('detection artifact');


% extract value of sigma for initial detection artifact
sigArt = estAll(3);
% distance at which the artifact exponential decreases past 0.1% of its
% amplitude at 1 frame
d01percent = -sigArt*log(0.001)+1;
disp(['0.1% distance: ',num2str(d01percent), 'frames']);



%% ========================================================================
%                   define following data cutoffs
% ========================================================================

% initial cutoff (fast lifetimes that are cut off because they may
% represent noisy detection artifact) - in frames
cutoff_start = ceil(d01percent);
% NOTE: I allow the option of using DIFFERENT values for cutoff_start in
% the fast and slow movies further below

% end cutoff (long lifetimes close to the length of the movie have noisy
% frequencies because of bad statistics) - in percent, e.g. 0.1 for 10%
cutoff_slow = 0.1;



%% ========================================================================
% 
%       determine which entries in the data structure represent fast or
%       slow movies and store the positions
%
% =========================================================================

% threshold for fast versus slow data (in seconds)
speed_thresh = 1;

% default values for later interpolation (in seconds)
ipFramerate_fast = 0.4;
ipFramerate_slow = 2;

% initialize variables
ct_fast = 1;
ct_slow = 1;
mintimespan_fast = 10000;
mintimespan_slow = 10000;

% In this implementation, since we want to average over the data without
% any breaks due to averaging over different numbers of histograms, the
% histograms are all cut to the same length. Thus, we need to know what's 
% the maximum common time span for all the fast/slow movies - as we only
% input movies with the same framerate (2s slow, 0.4s fast), the knowledge
% of the last frame suffices

% loop over all entries in the data structure
for i=1:length(data)
    % read current detection frequency
    framerate_i = data(i).framerate;
    % read current movie timespan - the framespan will later be reduced by
    % cutoffs, so it has to be modified
    framenum = data(i).movieLength;
    timespan = framerate_i*framenum;
    %framespan = floor(framenum*(1-cutoff_end)) - cutoff_start;
    %timespan = framerate_i*framespan;
    
    
    
    %  if movie is fast and data exists
    if ( (framerate_i<=speed_thresh) & ~isempty(data(i).lftHist) )
        % enter i as a position of a fast movie
        positions_fast(ct_fast) = i;
        ct_fast = ct_fast+1;
        mintimespan_fast = min(mintimespan_fast,timespan);
    %  if movie is slow and data exists    
    elseif ( (framerate_i>speed_thresh) & ~isempty(data(i).lftHist) )
        % enter i as a position of a slow movie
        positions_slow(ct_slow) = i;
        ct_slow = ct_slow+1;
        mintimespan_slow = min(mintimespan_slow,timespan);
         
    end % of if
end

% NOTE: now positions_fast contains the positions of fast movies,
% positions_slow the positions of the slow movies in data
        

% normalization window for later data (in seconds)
normStart_fast = 6;
normEnd_fast = (5*ipFramerate_fast)*floor(mintimespan_fast/(5*ipFramerate_fast));
normStart_slow = 30;
normEnd_slow = min(300,(5*ipFramerate_slow)*floor(mintimespan_slow/(5*ipFramerate_slow)));
% normalization window for later data (in frames)
normStartp_fast = round(normStart_fast/ipFramerate_fast);
normEndp_fast   = round(normEnd_fast/ipFramerate_fast);
normStartp_slow = round(normStart_slow/ipFramerate_slow);
normEndp_slow   = round(normEnd_slow/ipFramerate_slow);



%% ========================================================================
% 
%       Preparatory STAGE: average all fast/slow data 
%
% =========================================================================


% the following procedure is performed twice, once for fast, once for slow 
% movies; the respective interpolation time vectors were defined above


for r=1:2
    
    if r==1
        usepos      = positions_fast;
        ipFramerate = ipFramerate_fast;
        normw1      = normStartp_fast;
        normw2      = normEndp_fast;
        xlimit      = 60;
        cutoff_start= cutoff_start;
        maxl        = floor(mintimespan_fast/ipFramerate_fast);
        cutoff_end  = 0.1;
    elseif r==2
        usepos      = positions_slow;
        ipFramerate = ipFramerate_slow;
        normw1      = normStartp_slow;
        normw2      = normEndp_slow;
        xlimit      = 300;
        cutoff_start= cutoff_start;
        maxl        = floor(mintimespan_slow/ipFramerate_slow);
        cutoff_end  = cutoff_slow;
    end

    lup = length(usepos);    
    histMatrix  = nan*zeros(lup,1000);
    numPersistVec = nan*zeros(lup,1);
    
    figure
    hold on
    
    % loop over the appropriate positions in data
    for p=1:length(usepos)
        
        % index position of this movie
        i = usepos(p);
        % read current detection frequency
        currFramerate = data(i).framerate;  
        % read lifetime Histogram
        currHist = data(i).lftHist;
        % corresponding time vector
        tvec_curr = currFramerate*[1:length(currHist)];
         
        % determine number of cells/trajectories in this movie
        ncells(i) = sum(currHist(cutoff_start:length(currHist)));
        
        % correct histogram for movie length
        currCorrVec=lftHist_correctionVector(length(currHist));
        currLen = min( length(currHist),length(currCorrVec) );
        currHistCorr = currHist(1:currLen).*currCorrVec(1:currLen);
    
        
        % cut off points at the very end, because of bad counting
        % statistics at the end of the movie; the percentage of points at
        % the end to be cut off is determined above by the varaiable
        % cutoff_end
        cl = length(currHistCorr);
        cutpoint = round(maxl*(1-cutoff_end));
        currHistCorr( cutpoint:cl ) = [];
        tvec_curr( cutpoint:length(tvec_curr) ) = [];
        
        % perform cutoff at beginning, i.e. reject the first points because 
        % of the tracking/detection artifact
        currHistCorr(1:cutoff_start) = [];
        tvec_curr(1:cutoff_start) = [];
        
        % for 'slow' framerate, read number of persistent objects
        if (r==2) 
            % number of persistent trajectories = those longer than
            % cutpoint; if the lifetime vector/histogram of the cutoff
            % trajectories alreaday exists, use it, otherwise determine the
            % number separately
            csep = 0;
            if isfield(data,'lftHistCut')
                cuthist = data(i).lftHistCut;
                if ~isempty(cuthist)
                    numPers = sum( cuthist(cutpoint:length(cuthist)) );
                else
                    csep = 1;
                end
            else
                csep = 1;
            end
            
            if csep==1
                % if it isn't available yet, determine with function, but
                % note that since the threshold is in seconds, it's
                % necessary to multiply the cutoff point by framerate
                [numPers]=numPersistentField(data(i), currFramerate*cutpoint);
            end
            
            % normalize by the total sum, where the sum of the histogram plus
            % the number of persistent objects equals 100%
            totalsum = sum(currHistCorr)+numPers;
            ncells(i) = ncells(i) + numPers;
        else
            p120 = find(tvec_curr==120);
            totalsum = sum(currHistCorr(1:p120));
            
            totalsum = sum(currHistCorr);
            
        end
               
        % pre-normalize first by remaining sum
        currHistNorm = currHistCorr/totalsum;
                   
        % enter interpolated/renormalized histogram into summarizing matrix
        histMatrix(p,:) = nan;
        histMatrix(p,1:length(currHistNorm)) = currHistNorm;
        
                
        % plot results
        plot(tvec_curr,cumsum(currHistNorm),'b.-');
        
        if r==1, axis([0 120 0 1.05]); end
        if r==2, axis([0 600 0 1.05]); end
        
        %text( (0.2*xlimit), (0.9*max(currHistNorm)), ['movie # ',num2str(i)] );
        pause(0.3); 
                
        
    end % of for p
    
    % when the averaging takes place, that's a possible source for
    % artifacts for the first rayleigh population -
    % e.g. when the first defined point is averaged over 1 entry
    % and then the next one over 6 points - the first point is only if
    % the number of entries is >1 for more than 2 rows
    for ci = 1:20
        col = histMatrix(:,ci);
        n_nan = sum(isnan(col));
        n_def = lup - n_nan;
        if (n_nan>n_def) & (n_def>0)
            histMatrix(:,ci) = nan;
        end
    end
                   
              
    
    % calculate average
    histVectorAve = nanmean(histMatrix,1);
            
    % now find all finite point positions
    finitepos = find(isfinite(histVectorAve));
    
    % define appropriate restricted time vector
    tvecComp = tvec_curr(1)+ipFramerate*[1:length(histVectorAve)];
    tvecRes = tvecComp(finitepos);
    % define appropriate restricted lft vector
    histVectorRes = histVectorAve(finitepos);
    
    
    % if necessary crop more points here
    cp = 0;
    if cp>0, histVectorRes(1:cp)  = nan; end
    
%% continue HERE
    % if there's more than one entry for this condition
    if p>1
        histVectorResDel = [];
        for c=1:p
            histMatrixDel = histMatrix;
            % delete this rwo from the data matrix
            histMatrixDel(c,:) = [];
            % average
            histVectorAveDel = nanmean(histMatrixDel,1);
            % restrict to same positions as above
            histVectorResDel(c,:) = histVectorAveDel(finitepos);
        end % of for c
    end % of if
    
    
    if r==1
        histVector_fast = [ round(tvecRes*10)/10; histVectorRes ; histVectorResDel];
        totalNumCells_fast = sum(ncells(usepos));
        results.numcells_fast = totalNumCells_fast;
        plot(tvecRes, cumsum(histVectorRes),'r.-');
        
        figure;
        plot(tvecRes, histVectorRes,'r.-');
        axis([0 120 0 1.05*max(histVectorRes)]);
        
    elseif r==2
        histVector_slow = [ round(tvecRes*10)/10; histVectorRes ; histVectorResDel]; 
        totalNumCells_slow = sum(ncells(usepos));
        results.numcells_slow = totalNumCells_slow;
        plot(tvecRes, cumsum(histVectorRes),'r.-');
        
        figure;
        plot(tvecRes, 1-cumsum(histVectorRes),'r.-');
        guessvec = [0 0.3 4 2 0.3 8 2 0.3 90 1];
        fixvec = [0 0 0 0 0 0 0 0 0 0];
        [estSlowInv] = fitcurveMultiWeibullODF_lsq(tvecRes, 1-cumsum(histVectorRes), guessvec, fixvec);
        disp(['immobile fraction offset = ',num2str(estSlowInv(1))]);
              
        figure;
        plot(tvecRes, histVectorRes,'r.-');
        axis([0 600 0 1.05*max(histVectorRes)]);
        
        
    end % of if
    
    
end



%% write out results
% lftHist_fast = histVector_fast;
% lftHist_slow = histVector_slow;

results.hist_fast = histVector_fast;
results.hist_slow = histVector_slow;
% results.fitpar_fast = estFast;
% results.fitpar_slow = estSlow2;


end % of function




%% =============== mini subfunction for box filter
function [filteredtrace]=sboxfilter_edgefix(trace,boxlength)
len = length(trace);
for i=1:len
    cbox = min([(i-1),(len-i),boxlength]);
    tbeg = max(1,i-cbox);
    tend = min(len,i+cbox);
    filteredtrace(i)=mean(trace(tbeg:tend));
end
end % of subfunction



    
