function plusTipSettingSweepPlots(allData,saveDir)
% plusTipSettingSweepPlots displays allData in four plots per control parameter

if nargin<2 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Choose directory for storing allData plots');
end

% change name of maxTgapFrames to timeWindow so it matches first column
allData{1,cellfun(@(x) isequal(x,'maxTgapFrames'),allData(1,:))}='timeWindow';

% get unique parameter names that were tested from first column
paramsTested=unique(allData(2:end,1));


for i=1:length(paramsTested)
    % row indices for the current parameter
    r=find(cellfun(@(x) strcmpi(x,paramsTested{i}),allData(:,1)));
    paramInfo.(paramsTested{i}).rows=r;

    % column where the values tested for current parameter
    c=find(cellfun(@(x) strcmpi(x,paramsTested{i}),allData(1,:)));
    % here are the actual values
    paramInfo.(paramsTested{i}).xvals=cell2mat(allData(r,c));
end


statNames=allData(1,:)';
% remove all columns containing SEM values
idx2rem=find(cellfun(@(x) ~isempty(strfind(x,'_2')),statNames));
% also remove the first 9 columns (info about param values used)
idx2rem=[1:9 idx2rem'];
allData(:,idx2rem)=[];
% here are the stat names now available
statNames=allData(1,:);
% pretty up the names
statNames=cellfun(@(x) strrep(x,'_SE',''),statNames,'uniformoutput',0);
statNames=cellfun(@(x) strrep(x,'_','-'),statNames,'uniformoutput',0);
allData(1,:)=statNames;




% indices for for allData columns where 'percent' is in the name
prcntStatIdx=find(cellfun(@(x) ~isempty(strfind(x,'percent')),statNames));
% indices for growth/fgap/bgap speed mean
speedMean=find(cellfun(@(x) ~isempty(strfind(x,'speed-mean')),statNames));
% indices for growth/fgap/bgap lfetime mean
lifeMean=find(cellfun(@(x) ~isempty(strfind(x,'lifetime-mean')),statNames));
% indices for growth/fgap/bgap length mean (but not freq)
lengthMean=find(cellfun(@(x) ~isempty(strfind(x,'length-mean')),statNames) & cellfun(@(x) isempty(strfind(x,'freq')),statNames));
% indices for numbers of growths/fgaps/bgaps
numEvents=[find(cellfun(@(x) ~isempty(strfind(x,'nGrowths')) | ~isempty(strfind(x,'nFgaps')) | ~isempty(strfind(x,'nBgaps')),statNames)) find(cellfun(@(x) ~isempty(strfind(x,'nTracks')),statNames))];

% iterate thru control parameters
for i=1:length(paramsTested);

    r=paramInfo.(paramsTested{i}).rows;
    x=paramInfo.(paramsTested{i}).xvals;

    y={speedMean; lifeMean; lengthMean; numEvents; prcntStatIdx};
    
    for iPlot=1:length(y)
        mkPlot(allData,r,x,y{iPlot},paramsTested{i});
        
        switch iPlot
            case 1
                str='speedMean';
            case 2
                str='lifeMean';
            case 3
                str='lengthMean';
            case 4
                str='numEvents';
            case 5
                str='prcntStatIdx';
        end
        saveas(gcf,[saveDir filesep paramsTested{i} '_' str '.fig'])
        saveas(gcf,[saveDir filesep paramsTested{i} '_' str '.png'])
        close(gcf)

    end
end




% make plot
function mkPlot(allData,r,x,y,paramName)
% r=rows corresponding to current parameter
% x=column with values tested for the current parameter (x-axis values)
% y=columns of allData where data should be pulled for plotting

% get just the labels we need for current plot
ylabels=allData(1,y);

y=cell2mat(allData(r,y)); % y-axis values
xrep=repmat(x,[1 length(ylabels)]); % x-axis values

% properties for plotting
nLines=size(y,2);
nRep=ceil(nLines/3);
mrkTpe=repmat(['o';'^';'s'],nRep,1);
prop_name(1) = {'Marker'};
prop_name(2) = {'MarkerFaceColor'};
prop_name(3) = {'MarkerEdgeColor'};
prop_name(4) = {'Color'};

% set the colormap
cM=varycolor(nLines);
cMap=mat2cell(cM,ones(nLines,1),3);

figure; hold on
for iLine=1:nLines
    prop_values(1,1) = {mrkTpe(iLine)};
    prop_values(1,2) = cMap(iLine,:);
    prop_values(1,3) = cMap(iLine,:);
    prop_values(1,4) = cMap(iLine,:);

    xCoord=xrep(:,iLine);
    yCoord=y(:,iLine);

    h=plot(xCoord,yCoord);
    set(h,prop_name,prop_values)

    clear prop_values

end

legend(ylabels,'location','bestoutside')
xlabel(paramName)
set(gca,'XTick',x)
xlim([min(x) max(x)])


% close all