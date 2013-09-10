%postprocessing of tracksNA from ...
% [tracksNA,forceFC,forceFA,forceBGband,forceBG] = colocalizationAdhesionsWithTFM('/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/130429_cell11_100f/ROI','Region1',17,700,true);

% plot between force in NA on their birth vs. force in the BG
clearvars -except tracksNA forceFC forceFA forceBGband forceBG
nFrames  = tracksNA(1).iFrame(end);
forceEmergingNA = [];
forceBGbandNA = [];
forceBGNA = [];
p=0;
idx = false(numel(tracksNA),1);
for k=1:numel(tracksNA)
    % look for tracks that had a state of 'BA' and become 'NA'
    firstNAidx = find(strcmp(tracksNA(k).state,'NA'),1,'first');
    % see if the state is 'BA' before 'NA' state
    if (~isempty(firstNAidx) && firstNAidx>1 && strcmp(tracksNA(k).state(firstNAidx-1),'BA')) || (~isempty(firstNAidx) &&firstNAidx==1)
        p=p+1;
        idx(k) = true;
        tracksNA(k).emerging = true;
        tracksNA(k).emergingFrame = firstNAidx;
        forceEmergingNA(p) = tracksNA(k).forceMag(firstNAidx);
        forceBGbandNA(p) = forceBGband(firstNAidx).mean;
        forceBGNA(p) = forceBG(firstNAidx).mean;
    else
        tracksNA(k).emerging = false;
    end        
end
% [h,p,ci,stats] = ttest(forceEmergingNA,forceBGbandNA)
% h =  1
% p =   7.3323e-11
% ci =   45.5826   81.5166
% stats = 
%     tstat: 6.9856
%        df: 159
%        sd: 115.0716
%% Plotting
figure, notBoxPlot(forceBGNA,1)
notBoxPlot(forceBGbandNA,2)
notBoxPlot(forceEmergingNA,3)
%% Analysis of those whose force was under noise level: how long does it take
% Analysis shows that force is already developed somewhat compared to
% background. 
% For those NAs, first just see their force profile -> they are mostly
% increasing!
figure, hold on
meanNAForce = mean(forceEmergingNA);
% Filter out any tracks that has 'Out_of_ROI' in their status (especially
% after NA ...)
trNAonly = tracksNA(idx);
p=0; q=0;
xlim([0 60])
ylim([0 700])

for k=1:numel(tracksNA)
    if tracksNA(k).emerging && tracksNA(k).forceMag(tracksNA(k).emergingFrame)<2*meanNAForce ...
            && ~any(strcmp(tracksNA(k).state(tracksNA(k).emergingFrame:end),'Out_of_ROI'))
        % maturing NAs
        if (any(strcmp(tracksNA(k).state(tracksNA(k).emergingFrame:end),'FC')) || ...
                any(strcmp(tracksNA(k).state(tracksNA(k).emergingFrame:end),'FA'))) && ...
                sum(tracksNA(k).presence)>8
            
            tracksNA(k).maturing = true;
            p=p+1;
            t = 5/60*(tracksNA(k).iFrame(tracksNA(k).presence));
            f = tracksNA(k).forceMag(tracksNA(k).presence);
            [slopeMaturingR(p),slopeMaturing(p)] = regression(t,f);
            if slopeMaturing(p)>1.5
                plot(tracksNA(k).iFrame-tracksNA(k).emergingFrame,tracksNA(k).forceMag,'b','Linewidth',2)
            else
                plot(tracksNA(k).iFrame-tracksNA(k).emergingFrame,tracksNA(k).forceMag,'r','Linewidth',2)
                display(num2str(k))
            end
            forceMaturing{p} = tracksNA(k).forceMag(tracksNA(k).presence);
        elseif sum(tracksNA(k).presence)<61 && sum(tracksNA(k).presence)>6
        % failing NAs
            tracksNA(k).maturing = false;
            q=q+1;
            [slopeFailingR(q),slopeFailing(q)] = regression(5/60*(tracksNA(k).iFrame(tracksNA(k).presence)),tracksNA(k).forceMag(tracksNA(k).presence));
            if slopeFailing(q)>2.5
                plot(tracksNA(k).iFrame-tracksNA(k).emergingFrame,tracksNA(k).forceMag,'Linewidth',2,'Color',[0 .2 .2])
            elseif slopeFailing(q)>0
                plot(tracksNA(k).iFrame-tracksNA(k).emergingFrame,tracksNA(k).forceMag,'Linewidth',1,'Color',[.2 .2 .2])
            else
                plot(tracksNA(k).iFrame-tracksNA(k).emergingFrame,tracksNA(k).forceMag,'Linewidth',2,'Color',[.5 .1 .1])
%                 display(num2str(k))
            end
            forceFailing{q} = tracksNA(k).forceMag(tracksNA(k).presence);
        end
    end
end
%% Plotting
figure, notBoxPlot(slopeMaturing,1)
notBoxPlot(slopeFailing,2)
slopeMaturing=slopeMaturing';
slopeFailing=slopeFailing';
slopeMaturingMean = mean(slopeMaturing);
slopeMaturingErr = std(slopeMaturing)/sqrt(length(slopeMaturing));
slopeFailingMean = mean(slopeFailing);
slopeFailingErr = std(slopeFailing)/sqrt(length(slopeFailing));
%% average profile for the first n frames
nSampleFrames = 25;
forceMaturingArray = NaN(numel(forceMaturing),nSampleFrames);
forceFailingArray = NaN(numel(forceFailing),nSampleFrames);
for k=1:numel(forceMaturing)
    curProfileMaturing = forceMaturing{k};
    fmax = min(nSampleFrames, length(curProfileMaturing));
    forceMaturingArray(k,1:fmax) = curProfileMaturing(1:fmax);
end
for k=1:numel(forceFailing)
    curProfileFailing = forceFailing{k};
    fmax = min(nSampleFrames, length(curProfileFailing));
    forceFailingArray(k,1:fmax) = curProfileFailing(1:fmax);
end
maturingFavg = nanmean(forceMaturingArray,1)';
failingFavg = nanmean(forceFailingArray,1)';
maturingFstd = nanstd(forceMaturingArray,1)';
failingFstd = nanstd(forceFailingArray,1)';
forceFailingArray = forceFailingArray';
forceMaturingArray = forceMaturingArray';
%% Analysis of NAs that becomes FAs with their area and force
