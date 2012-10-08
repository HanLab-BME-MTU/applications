function plusTipQuadScatter(xAxisInfo,yAxisInfo,groupList,remBegEnd,timeRange,doPlot,saveDir)
% Makes quadrant scatterplot split by property percentiles or values, mapped onto image

% SYNOPSIS: plusTipQuadScatter(xAxisInfo,yAxisInfo,groupList,remBegEnd,timeRange,doPlot,saveDir)

% INPUT:    xAxisInfo/yAxisInfo  : structures containing the following fields:
%               .name
%                   MT dynamics parameter, chosen from the following list,
%                   where both are from the same group (i.e. both from growth):
%                            'growthSpeed' (track average speed, microns/min)
%                            'growthLifetime' (growth lifetime, frames)
%                            'growthDisp' (displacement during growth, microns)
%                            'fgapSpeed' (etc...)
%                            'fgapLifetime'
%                            'fgapDisp'
%                            'bgapSpeed'
%                            'bgapLifetime'
%                            'bgapDisp'
%               .splitPercentile
%                   percentile where the data should be split for the
%                   particular parameter
%               .splitValue
%                   values where the data should be split for the
%                   particular parameter.
%               .minMax
%                   [lower upper] limits to x and y axes of scatter plot
%                   data outside this range will be discarded.
%
%                NOTE: only percentile OR splitValue value should be
%                      given for each parameter - the other will be
%                      calculated accordingly
%
%           groupList: output of plusTipPickGroups, containing groups of
%                      projects to be analyzed together
%           remBegEnd      : 1 to remove tracks existing at the beginning
%                            or end of the timeRange
%           timeRange      : frame range over which data should be used
%           doPlot         : (1) to make all plots, 2 to only make summary
%                            plots
%
% OUTPUT:
%           A scatterplot of the parameters plotted against each other,
%           split into four colors at the boundaries created by the
%           division marks, and an image overlay where the colors
%           correspond to the tracks of interest.  Four individual overlays
%           are also created, displaying each color in turn.  A percentage
%           bar showing the relative proportion of the subpopulations is
%           saved as a separate figure.
%
%           If the same property is given for both parameters,
%           there will be 2-3 colors splitting the data into only 2 or 3
%           groups.  This can be done to show the first, fourth, and middle
%           two quartiles of, say, growth speed, if cutoffs 1 and 2 are 25
%           and 75, respectively.
%
%           If there are mutliple groups, a stacked percentage bar will
%           also be saved.  btwGrpQuadStats.mat is a structure containing
%           the parameters used as well as the track counts/percentages in
%           each subpopulation for each movie.
%
% 2009/08 - Kathryn Applegate Matlab R2008a

if nargin<2 || ~isstruct(xAxisInfo) || ~isstruct(yAxisInfo)
    msgbox('plusTipQuadScatter: not enough input parameters')
    return
end

xFields=fieldnames(xAxisInfo);
yFields=fieldnames(yAxisInfo);
if ~isequal(sort(xFields),sort({'name';'splitPercentile';'splitValue';'minMax'})) ||...
        ~isequal(sort(yFields),sort({'name';'splitPercentile';'splitValue';'minMax'}))
    msgbox('plusTipQuadScatter: one or more fieldnames for xAxisInfo or yAxisInfo are not correct')
    return
end

if ~isempty(xAxisInfo.minMax)
    %xAxisInfo.minMax=xAxisInfo.minMax;
else
    xAxisInfo.minMax=[-inf inf];
end
if ~isempty(yAxisInfo.minMax)
    %yAxisInfo.minMax=yAxisInfo.minMax;
else
    yAxisInfo.minMax=[-inf inf];
end

if ~xor(isempty(xAxisInfo.splitValue),isempty(xAxisInfo.splitPercentile)) ||...
        ~xor(isempty(yAxisInfo.splitValue),isempty(yAxisInfo.splitPercentile))
    msgbox('plusTipQuadScatter: splitValue or splitPercentile value should be empty','Input Error','error')
    return
end

% check projData input
if nargin<3 || isempty(groupList)
    [groupList]=combineGroupListFiles(0);
end

% assume we should plot all data
if nargin<4 || isempty(remBegEnd)
    remBegEnd=0;
end

% check whether a time range for plotting was input
if nargin<5 || isempty(timeRange)
    timeRange=[1 inf];
end

if nargin<6 || isempty(doPlot)
    doPlot=1;
end

if nargin<7 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Please choose output directory');
end
if ~isdir(saveDir)
    mkdir(saveDir)
end

% define which parameter sets are complemetary.  xAxisInfo.name and
% yAxisInfo.name must come from the same group
g1 = {'growthSpeed', 'growthLifetime', 'growthDisp'};
g2 = {'fgapSpeed', 'fgapLifetime', 'fgapDisp'};
g3 = {'bgapSpeed', 'bgapLifetime', 'bgapDisp'};

% check to make sure they are from the same group
errFlag=1;
if ismember(xAxisInfo.name,g1) && ismember(yAxisInfo.name,g1)
    trackType=1;
    errFlag=0;
end
if ismember(xAxisInfo.name,g2) && ismember(yAxisInfo.name,g2)
    trackType=2;
    errFlag=0;
end
if ismember(xAxisInfo.name,g3) && ismember(yAxisInfo.name,g3)
    trackType=3;
    errFlag=0;
end
if errFlag==1
    msgbox('plusTipQuadScatter: Input parameters must be from same group (i.e. growth)','Input Error','error')
    return
end

fileExt='.tif';


% get the group names
isML = isa(groupList,'MovieList');
if isML
    [~,grpNames] = arrayfun(@(x) fileparts(x.getFilename()), groupList,...
        'Unif', false);
else
    [grpNames,b,m]=unique(groupList(:,1));
    [b,idx]=sort(b);
    grpNames=grpNames(idx);
end

btwGrpQuadStats.grpPopRYGB=zeros(length(grpNames),4);
btwGrpQuadStats.grpPrctRYGB=zeros(length(grpNames),4);

% iterate through each group
nProj=size(groupList,1);
nGrps=length(grpNames);

s = length(num2str(nProj));
iPrjStr = sprintf('%%.%dd',s);
s = length(num2str(nGrps));
iGrpStr = sprintf('%%.%dd',s);

if nProj>1
    pVSpDir=[saveDir filesep xAxisInfo.name 'VS' yAxisInfo.name];
else
    pVSpDir=[saveDir filesep 'quad' filesep xAxisInfo.name 'VS' yAxisInfo.name];
end
if isdir(pVSpDir)
    rmdir(pVSpDir,'s')
end

data1all=[];
data2all=[];
for iGroup=1:nGrps

    groupNum = sprintf(iGrpStr,iGroup);
    grpDir=[pVSpDir filesep 'grp' groupNum '_' grpNames{iGroup}];
    mkdir(grpDir)

    if ~isML
        % indices of projects in iGroup
        grpIdx=find(cellfun(@(x) ~isempty(strmatch(x,grpNames{iGroup},'exact')),groupList(:,1)));
        nProj =length(grpIdx);
    else
        nProj = numel(groupList(iGroup).getMovies());
    end
    

    % initialize percentile, threshold, and population arrays
    quadStats.splitPrctlX=zeros(nProj, 1);
    quadStats.splitValueX=zeros(nProj, 1);
    quadStats.splitPrctlY=zeros(nProj, 1);
    quadStats.splitValueY=zeros(nProj, 1);
    quadStats.grpPopRYGB=zeros(nProj, 4);
    quadStats.grpPrctRYGB=zeros(nProj, 4);

    data1grp=[];
    data2grp=[];
    % load data and make plots for each project iSub in group iGroup
    for iSub=1:nProj

        % load data for iSub
        if ~isML
            projDir=groupList{grpIdx(iSub),2};
            load([projDir filesep 'meta' filesep 'projData'])
            
            % find actual frame range for this movie (smaller if range input)
            tempTimeRange=timeRange;
            tempTimeRange(1)=max([tempTimeRange(1);projData.detectionFrameRange(1);...
                projData.trackingFrameRange(1);projData.postTrackFrameRange(1)]);
            tempTimeRange(2)=min([tempTimeRange(2);projData.detectionFrameRange(2);...
                projData.trackingFrameRange(2);projData.postTrackFrameRange(2)]);
            
            % get first image from imDir
            [listOfImages] = searchFiles('.tif',[],formatPath(projData.imDir),0);
            img = double(imread([listOfImages{1,2} filesep listOfImages{1,1}]));
        else
            movie = groupList(iGroup).getMovies{iSub};
            
            iDetProc = movie.getProcessIndex('DetectionTrackingProcess',1,0);
            detProc = movie.getProcess(iDetProc);
            if isa(detProc, 'CometDetectionProcess');
                detectionRange = [detProc.funParams_.firstFrame detProc.funParams_.lastFrame];
            else
                detectionRange = [1 movie.nFrames_];
            end
            tempTimeRange(1)=max([timeRange(1); detectionRange(1)]);
            tempTimeRange(2)=min([timeRange(2); detectionRange(2)]);
            
            % Read post-tracking info
            iPostProc = movie.getProcessIndex('CometPostTrackingProcess',1,0);
            postProc = movie.getProcess(iPostProc);
            iChan = find(postProc.checkChannelOutput,1);
            projData= postProc.loadChannelOutput(iChan,'output','projData');
            
            % Read first image
            img = movie.channels_(iChan).loadImage(1);
        end




        % get size of the figures
        minY=1; maxY=size(img,1);
        minX=1; maxX=size(img,2);
        scrsz = get(0,'ScreenSize');
        screenW=scrsz(3);
        screenL=scrsz(4);
        magCoef=inf;
        maxMagCoefW = (0.8*screenW)/(maxX-minX+1);
        maxMagCoefL = (0.8*screenL)/(maxY-minY+1);
        if magCoef > min([maxMagCoefW; maxMagCoefL])
            calcMagCoef = min([magCoef; maxMagCoefW; maxMagCoefL]);
        else
            calcMagCoef = magCoef;
        end
        movieL = (calcMagCoef*(maxY-minY+1));
        movieW = (calcMagCoef*(maxX-minX+1));
        figPos=[round(screenW*(1-movieW/screenW)/2) round(screenL*(1-movieL/screenL)/2) movieW movieL];


        % extract merged track info (block matrix) - take abs value and
        % convert displacement to microns and lifetimes to sec
        dataMatMerge=plusTipMergeSubtracks(projData); % merged data is first output
        allData=abs(dataMatMerge);
        allData(:,6)=allData(:,6).*projData.secPerFrame;
        allData(:,7)=allData(:,7).*(projData.pixSizeNm./1000);
        % get all the subtrack coordinates corresponding to the merged data
        [xMat,yMat]=plusTipGetSubtrackCoords(projData,[],1);

        % for xAxisInfo.name and yAxisInfo.name, get data and label
        for iParam=1:2
            % assign data based on user input
            if iParam==1
                paramName=xAxisInfo.name;
            else
                paramName=yAxisInfo.name;
            end

            switch paramName
                case 'growthSpeed'
                    tempData=allData(:,4);
                    tempLabel='Growth Speed (microns/min)';
                case 'fgapSpeed'
                    tempData=allData(:,4);
                    tempLabel='Fgap Speed (microns/min)';
                case 'bgapSpeed'
                    tempData=allData(:,4);
                    tempLabel='Bgap Speed (microns/min)';
                case 'growthLifetime'
                    tempData=allData(:,6);
                    tempLabel='Growth Lifetime (sec)';
                case 'fgapLifetime'
                    tempData=allData(:,6);
                    tempLabel='Fgap Lifetime (sec)';
                case 'bgapLifetime'
                    tempData=allData(:,6);
                    tempLabel='Bgap Lifetime (sec)';
                case 'growthDisp'
                    tempData=allData(:,7);
                    tempLabel='Growth Displacement (microns)';
                case 'fgapDisp'
                    tempData=allData(:,7);
                    tempLabel='Fgap Displacement (microns)';
                case 'bgapDisp'
                    tempData=allData(:,7);
                    tempLabel='Bgap Displacement (microns)';
            end

            if iParam==1
                data1=tempData;
                label1=tempLabel;
                tempLabel1=tempLabel;
            else
                data2=tempData;
                label2=tempLabel;
                tempLabel2=tempLabel;
            end
        end

        if remBegEnd==1
            % any track not entirely contained within the frame range will be excluded
            subtrackIdx=(allData(:,5)==trackType & allData(:,2)>tempTimeRange(1) & allData(:,3)<tempTimeRange(2));
        else
            % any track which ends before the frame range begins or begins after the frame range ends will be excluded.
            subtrackIdx=(allData(:,5)==trackType & allData(:,2)<tempTimeRange(2) & allData(:,3)>tempTimeRange(1));
        end
        inRangeIdx=data1>=xAxisInfo.minMax(1) & data1<=xAxisInfo.minMax(2) & data2>=yAxisInfo.minMax(1) & data2<=yAxisInfo.minMax(2);


        % remove out-of-range values from x-axis data
        data1=data1(subtrackIdx & inRangeIdx);
        if isempty(xAxisInfo.splitValue)
            splitPercentileX=xAxisInfo.splitPercentile;
            splitValueX=prctile(data1,splitPercentileX);
        else
            splitValueX=xAxisInfo.splitValue;
            splitPercentileX=floor(100*length(find(data1<splitValueX))/length(data1));
        end
        label1=[label1 ', P' num2str(splitPercentileX) ' = ' sprintf('%3.2f',splitValueX)];
        % remove out-of-range values from y-axis data
        data2=data2(subtrackIdx & inRangeIdx);
        if isempty(yAxisInfo.splitValue)
            splitPercentileY=yAxisInfo.splitPercentile;
            splitValueY=prctile(data2,splitPercentileY);
        else
            splitValueY=yAxisInfo.splitValue;
            splitPercentileY=floor(100*length(find(data2<splitValueY))/length(data2));
        end
        label2=[label2 ', P' num2str(splitPercentileY) ' = ' sprintf('%3.2f',splitValueY)];

        xMat=xMat(subtrackIdx & inRangeIdx,:);
        yMat=yMat(subtrackIdx & inRangeIdx,:);


        % store all data in iGroup list
        data1grp = [data1grp; data1];
        data2grp = [data2grp; data2];
        % store data in all-group list
        data1all = [data1all; data1];
        data2all = [data2all; data2];

        colorMap=['b','g','y','r'];
        popRYGB=zeros(1,4);

        if doPlot==1
            h=zeros(8,1); % for figure handle storage
            % fig 6: make scatterplot
            h(6)=figure;

            % fig 1: overlay of all 4 colors
            h(1)=figure('Position',figPos);
            imagesc(img); colormap gray;
            hold on
        end
        for iColor=1:4
            switch iColor
                case 1 % blue
                    idx=find(data1>splitValueX & data2>splitValueY);
                    popRYGB(4)=length(idx);
                case 2 % green
                    idx=find(data1<=splitValueX & data2>splitValueY);
                    popRYGB(3)=length(idx);
                case 3 % yellow
                    idx=find(data1>splitValueX & data2<=splitValueY);
                    popRYGB(2)=length(idx);
                case 4 % red
                    idx=find(data1<=splitValueX & data2<=splitValueY);
                    popRYGB(1)=length(idx);
            end
            if doPlot==1
                % figs 2-5: b,g,y,r
                h(iColor+1)=figure('Position',figPos);
                imagesc(img); colormap gray;
                hold on
                plot(xMat(idx,:)',yMat(idx,:)',colorMap(iColor))
                axis equal

                figure(h(1)) % fig 1: overall all four colors
                plot(xMat(idx,:)',yMat(idx,:)',colorMap(iColor))
                axis equal

                figure(h(6)) % fig 6: scatterplot single
                scatter(data1(idx),data2(idx),[],colorMap(iColor),'.');
                hold on
            end
        end

        if doPlot==1
            % apply limits and add title to fig6
            figure(h(6))
            xlim(xAxisInfo.minMax)
            ylim(yAxisInfo.minMax)
            xlabel(label1)
            ylabel(label2)
            if remBegEnd==1
                % any track not entirely contained within the frame range will be excluded
                title({['N = ' num2str(length(data1)) ' tracks']; ...
                    ['starting after frame ' num2str(tempTimeRange(1))...
                    ' and ending before frame ' num2str(tempTimeRange(2))]});
            else
                % any track which ends before the frame range begins or begins after the frame range ends will be excluded.
                title({['N = ' num2str(length(data1)) ' tracks'];...
                    ['starting at or before frame ' num2str(tempTimeRange(2))...
                    ' and ending at or after frame ' num2str(tempTimeRange(1))]});
            end


            % make the percent bar for this movie
            [prctRYGB]=plusTipQuadColorbar(popRYGB);
            h(7)=gcf;
            xlabel(['Percents: red=' num2str(prctRYGB(1)) ', yellow=' ...
                num2str(prctRYGB(2)) ', green=' num2str(prctRYGB(3)) ', blue=' num2str(prctRYGB(4))])
            
            if ~isML
                projNum = sprintf(iPrjStr,grpIdx(iSub));
            else
                projNum = sprintf(iPrjStr,iSub);
            end
            tempStr=[grpDir filesep 'groupListIdx_' projNum];

            %saveas(h(6),[tempStr '_1_scatter' '.fig'])
            saveas(h(6),[tempStr '_1_scatter' fileExt])

            %saveas(h(7),[tempStr '_2_prctBar' '.fig'])
            saveas(h(7),[tempStr '_2_prctBar' fileExt])

            %saveas(h(1),[tempStr '_3_all' '.fig'])
            saveas(h(1),[tempStr '_3_all' fileExt])

            %saveas(h(2),[tempStr '_4_blue' '.fig'])
            saveas(h(2),[tempStr '_4_blue' fileExt])

            %saveas(h(3),[tempStr '_5_green' '.fig'])
            saveas(h(3),[tempStr '_5_green' fileExt])

            %saveas(h(4),[tempStr '_6_yellow' '.fig'])
            saveas(h(4),[tempStr '_6_yellow' fileExt])

            %saveas(h(5),[tempStr '_7_red' '.fig'])
            saveas(h(5),[tempStr '_7_red' fileExt])

            for i=1:7
                close(h(i))
            end
        else
            [prctRYGB]=plusTipQuadColorbar(popRYGB,doPlot);
        end

        quadStats.splitValueX(iSub,1)=splitValueX;
        quadStats.splitPrctlX(iSub,1)=splitPercentileX;
        quadStats.splitValueY(iSub,1)=splitValueY;
        quadStats.splitPrctlY(iSub,1)=splitPercentileY;
        quadStats.grpPopRYGB(iSub,:)=popRYGB;
        quadStats.grpPrctRYGB(iSub,:)=prctRYGB;
    end

    if nProj>1
        % make the percentage bar for all iGroup projects in the same plot
        plusTipQuadColorbar(quadStats.grpPopRYGB);
        if ~isML
            titleStr=strrep(groupList{grpIdx(iSub),1},'_','-');
        else
            titleStr=strrep([grpNames{iGroup} '_' num2str(iSub)],'_','-');
        end
            
        title(titleStr)
        %saveas(gcf,[pVSpDir '_prctBarAll' '.fig'])
        saveas(gcf,[grpDir filesep 'prctBarAll' fileExt])
        close(gcf)

        % sum the populations to make one percentage bar for the whole group
        merPop=sum(quadStats.grpPopRYGB,1);
        [merPrcl]=plusTipQuadColorbar(merPop);

        quadStats.grpPopMergedRYGB=merPop;
        quadStats.grpPrctMergedRYGB=merPrcl;

        btwGrpQuadStats.grpPopRYGB(iGroup,:)=merPop;

%         titleStr=strrep(groupList{grpIdx(iSub),1},'_','-');
        title([titleStr ', N=' num2str(nProj)])
        %saveas(gcf,[pVSpDir '_prctBarAll_Merged' '.fig'])
        saveas(gcf,[grpDir filesep 'prctBarAll_Merged' fileExt])
        close(gcf)


        % make the scatterplot for all data from the group
        h(8)=figure; % fig 8: scatterplot single
        for iColor=1:4
            switch iColor
                case 1 % blue
                    idx=find(data1grp>splitValueX & data2grp>splitValueY);
                case 2 % green
                    idx=find(data1grp<=splitValueX & data2grp>splitValueY);
                case 3 % yellow
                    idx=find(data1grp>splitValueX & data2grp<=splitValueY);
                case 4 % red
                    idx=find(data1grp<=splitValueX & data2grp<=splitValueY);
            end
            scatter(data1grp(idx),data2grp(idx),[],colorMap(iColor),'.');
            hold on
        end
        % apply limits and add title to fig8
        figure(h(8))
        xlim(xAxisInfo.minMax)
        ylim(yAxisInfo.minMax)
        xlabel(tempLabel1)
        ylabel(tempLabel2)
        if remBegEnd==1
            % any track not entirely contained within the frame range will be excluded
            title({['N = ' num2str(length(data1grp)) ' tracks']; ...
                ['starting after frame ' num2str(tempTimeRange(1))...
                ' and ending before frame ' num2str(tempTimeRange(2))]});
        else
            % any track which ends before the frame range begins or begins after the frame range ends will be excluded.
            title({['N = ' num2str(length(data1grp)) ' tracks'];...
                ['starting at or before frame ' num2str(tempTimeRange(2))...
                ' and ending at or after frame ' num2str(tempTimeRange(1))]});
        end
        saveas(h(8),[grpDir filesep 'scatterAll_Merged' fileExt])
        close(h(8))

        save([grpDir filesep 'allGrpQuadStats'],'quadStats','xAxisInfo','yAxisInfo');
    end
end

if nGrps>1
    [btwGrpQuadStats.grpPrctRYGB]=plusTipQuadColorbar(btwGrpQuadStats.grpPopRYGB);
    title(['All Groups, N=' num2str(nGrps)])
    %saveas(gcf,[pVSpDir filesep 'prctBarAllGroups' '.fig'])
    saveas(gcf,[pVSpDir filesep 'prctBarAllGroups' fileExt])
    close(gcf)
    
    
    % make the scatterplot for all data from all groups
        h(9)=figure; % fig 8: scatterplot single
        for iColor=1:4
            switch iColor
                case 1 % blue
                    idx=find(data1all>splitValueX & data2all>splitValueY);
                case 2 % green
                    idx=find(data1all<=splitValueX & data2all>splitValueY);
                case 3 % yellow
                    idx=find(data1all>splitValueX & data2all<=splitValueY);
                case 4 % red
                    idx=find(data1all<=splitValueX & data2all<=splitValueY);
            end
            scatter(data1all(idx),data2all(idx),[],colorMap(iColor),'.');
            hold on
        end
        % apply limits and add title to fig8
        figure(h(9))
        xlim(xAxisInfo.minMax)
        ylim(yAxisInfo.minMax)
        xlabel(tempLabel1)
        ylabel(tempLabel2)
        if remBegEnd==1
            % any track not entirely contained within the frame range will be excluded
            title({['N = ' num2str(length(data1all)) ' tracks']; ...
                ['starting after frame ' num2str(tempTimeRange(1))...
                ' and ending before frame ' num2str(tempTimeRange(2))]});
        else
            % any track which ends before the frame range begins or begins after the frame range ends will be excluded.
            title({['N = ' num2str(length(data1all)) ' tracks'];...
                ['starting at or before frame ' num2str(tempTimeRange(2))...
                ' and ending at or after frame ' num2str(tempTimeRange(1))]});
        end
        saveas(h(9),[pVSpDir filesep 'scatterAllGroups' fileExt])
        close(h(9))
    
    

    save([pVSpDir filesep 'btwGrpQuadStats'],'btwGrpQuadStats','xAxisInfo','yAxisInfo','groupList');
end

disp('Quadrant Scatter Plots...Finished!')


