function [thresh1,thresh2]=plusTipParamPlot(param1,cutoff1,param2,cutoff2,projData,remBegEnd)
% Makes quadrant scatterplot split by property percentiles, mapped  onto image

% SYNOPSIS: [thresh1,thresh2]=plusTipParamPlot(param1,cutoff1,param2,cutoff2,projData)

% INPUT:    param1/param2  : properties of MT dynamics, chosen from the
%                            following list, where both are from the same
%                            group (ie both from growth):
%                            'growthSpeed' (track average speed, microns/min)
%                            'growthLifetime' (growth lifetime, frames)
%                            'growthDisp' (displacement during growth, microns)
%                            'fgapSpeed'
%                            'fgapLifetime'
%                            'fgapDisp'
%                            'bgapSpeed'
%                            'bgapLifetime'
%                            'bgapDisp'
%           cutoff1/cutoff2: percentiles where the data should be split for
%                            properties param1/param2, respectively
%           projData (opt) : MT tracking data stored in the meta folder
%                            if not given, user will be asked for it
%           remBegEnd (opt): 1 to remove tracks existing at the beginning
%                            or end of the movie
%
% OUTPUT:   thresh1/thresh2: the parameter values which are the percentile
%                            values of the cutoffs
%           The principal outputs are two plots: one, a scatterplot of the
%           parameters plotted against each other, split into four colors
%           at the boundaries created by the percentiles, and the second,
%           an image overlay where the colors correspond to the tracks of
%           interest.  If the same property is given for both parameters,
%           there will be 2-3 colors splitting the data into only 2 or 3
%           groups.  This can be done to show the first, fourth, and middle
%           two quartiles of, say, growth speed, if cutoffs 1 and 2 are 25
%           and 75, respectively.
%
% 2009/08 - Kathryn Applegate Matlab R2008a


if nargin<4
    param1='growthSpeed';
    cutoff1=50;
    param2='growthLifetime';
    cutoff2=50;
end



% define which parameter sets are complemetary.  param1 and 2 must come
% from the same group
group1 = {'growthSpeed', 'growthLifetime', 'growthDisp'};
group2 = {'fgapSpeed', 'fgapLifetime', 'fgapDisp'};
group3 = {'bgapSpeed', 'bgapLifetime', 'bgapDisp'};

% check to make sure they are from the same group
errFlag=1;
if ismember(param1,group1) && ismember(param2,group1)
    groupNum=1;
    errFlag=0;
end
if ismember(param1,group2) && ismember(param2,group2)
    groupNum=2;
    errFlag=0;
end
if ismember(param1,group3) && ismember(param2,group3)
    groupNum=3;
    errFlag=0;
end
if errFlag==1
    errordlg('plusTipParamPlot: Input parameters must be from same group (i.e. growth)')
    uiwait
    return
end

% check projData input
if nargin<5 || isempty(projData)
    dirName=uigetdir(pwd,'Please select project directory');
    temp=load([dirName filesep 'meta' filesep 'projData.mat']);
    projData=temp.projData;
end

% assume we should plot all data
if nargin<6 || isempty(remBegEnd)
    remBegEnd=0;
end

% format image/analysis directory paths
imDir=formatPath(projData.imDir);
anDir=formatPath(projData.anDir);

% get first image from imDir
[listOfImages] = searchFiles('.tif',[],imDir,0);
img = double(imread([listOfImages{1,2} filesep listOfImages{1,1}]));

% extract track info (block matrix)
[dataMatMerge,dataMatReclass,percentFgapsReclass]=plusTipMergeSubtracks(projData);
allData=abs(dataMatMerge);

% get track type and matrices containing xy-coordinates for all subtracks
trackType=allData(:,5);
trackStarts=allData(:,2);
trackEnds=allData(:,3);
[xMat,yMat]=plusTipGetSubtrackCoords(projData,[],1);

if remBegEnd==1
    sF=min(allData(:,2));
    eF=max(allData(:,3));
else
    sF=0;
    eF=inf;
end

if ismember(groupNum,[1 2 3])
    subtrackIdx=find(trackType==groupNum & trackStarts>sF & trackEnds<eF);
end

% for param1 and param2, get data and label
for iParam=1:2
    % assign data based on user input    
    switch eval(['param' num2str(iParam)])

        case {'growthSpeed', 'fgapSpeed', 'bgapSpeed'} % speeds in microns/min
            % extract data
            tempData=allData(subtrackIdx,4);
            
            % assign axis labels
            switch eval(['param' num2str(iParam)])
                case 'growthSpeed'
                    tempLabel='Growth Speed (microns/min)';
                case 'fgapSpeed'
                    tempLabel='Fgap Speed (microns/min)';
                case 'bgapSpeed'
                    tempLabel='Bgap Speed (microns/min)';
            end

        case {'growthLifetime', 'fgapLifetime', 'bgapLifetime'} % lifetimes in seconds
            % extract data
            tempData=allData(subtrackIdx,6).*projData.secPerFrame;

            % assign axis labels
            switch eval(['param' num2str(iParam)])
                case 'growthLifetime'
                    tempLabel='Growth Lifetime (s)';
                case 'fgapLifetime'
                    tempLabel='Fgap Lifetime (s)';
                case 'bgapLifetime'
                    tempLabel='Bgap Lifetime (s)';
            end

        case {'growthDisp', 'fgapDisp', 'bgapDisp'} % displacements in microns
            % extract data
            tempData=(allData(subtrackIdx,7).*projData.pixSizeNm)./1000;

            % assign axis labels
            switch eval(['param' num2str(iParam)])
                case 'growthDisp'
                    tempLabel='Growth Displacement (microns)';
                case 'fgapDisp'
                    tempLabel='Fgap Displacement (microns)';
                case 'bgapDisp'
                    tempLabel='Bgap Displacement (microns)';
            end

        otherwise
            tempData='This parameter is not supported.';
    end

    % first iteration gives x-axis data, second gives y-axis data
    if iParam==1
        data1=tempData;
        thresh1=prctile(data1,cutoff1);
        label1=[tempLabel ', P' num2str(cutoff1) ' = ' sprintf('%3.2f',thresh1)];
    else
        data2=tempData;
        thresh2=prctile(data2,cutoff2);
        label2=[tempLabel ', P' num2str(cutoff2) ' = ' sprintf('%3.2f',thresh2)];
    end

end

colorMap=['b','g','y','r'];




% make scatterplot showing data distribution between the two parameters
figure 
hold on
for iColor=1:4
    switch iColor
        case 1 % blue
            idx=find(data1>thresh1 & data2>thresh2);
        case 2 % green
            idx=find(data1<=thresh1 & data2>thresh2);
        case 3 % yellow
            idx=find(data1>thresh1 & data2<=thresh2);
        case 4 % red
            idx=find(data1<=thresh1 & data2<=thresh2);     
    end

    scatter(data1(idx),data2(idx),[],colorMap(iColor),'.');
    xlabel(label1)
    ylabel(label2)
    title(['N = ' num2str(length(data1))]);
end

% plot the corresponding tracks on the image
h=figure; 
imagesc(img); colormap gray;
hold on
for iColor=1:4
    switch iColor
        case 1 % blue
            idx=find(data1>thresh1 & data2>thresh2);
        case 2 % green
            idx=find(data1<=thresh1 & data2>thresh2);
        case 3 % yellow
            idx=find(data1>thresh1 & data2<=thresh2);
        case 4 % red
            idx=find(data1<=thresh1 & data2<=thresh2);     
    end
    % make individual plot
    figure
    imagesc(img); colormap gray;
    hold on
    plot(xMat(subtrackIdx(idx),:)',yMat(subtrackIdx(idx),:)',colorMap(iColor))

    % make overlay
    figure(h)
    plot(xMat(subtrackIdx(idx),:)',yMat(subtrackIdx(idx),:)',colorMap(iColor))

end
