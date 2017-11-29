function [FCondMean,FCondPooled,FPerCell] = plotDifferentPerspectives( dataPerCell,grpIDsVect,varargin)
% plotDifferentPerspectives
% 
% dataPerCell: rxc matrix where r is the maximum number of measurements
% observed per cell and c is the number of cells (total in for all groups) 
% 
% grpIDsVect: rx1 vector designating the groupIDs for each of the columns
% in dataPerCell. 
% 

%% TEST : assume 100 dataPoints for each cell, 30 cells, 3 conditions, 10 cells per condition 
% dataPerCell = randn(100,30); 
% grpIDsCell = arrayfun(@(x) repmat(x,[10,1]),1:3,'uniformoutput',0); 
% grpIDsVect = vertcat(grpIDSCell{:}); 
% % note if you have a cell array for project paths per group you can
% generate the grpIDs by 
% grpIDsCell = cellfun(@(x) repmat(x,[length(x),1],1:numel(projList)); 
%% 
ip = inputParser;

ip.CaseSensitive = false;
ip.addParameter('names',[]);
ip.addParameter('colors',[]); 
ip.addParameter('colorShadesSEM',[]);
ip.addParameter('colorShadesSTD',[]); 
ip.parse(varargin{:});


%% Defaults for some nice colors
% List colors default
% black
% red
% organge
% pink
% blue
% nice cyan
% purple
% nice green 
%
if isempty(ip.Results.colors)
colors = {[0,0,0];[1,0,0];[0.996078431372549,0.698039215686275,0.298039215686275]; ...
    [0.905882352941177,0.160784313725490,0.541176470588235];[0.192156862745098,0.509803921568627,0.741176470588235]; ...
    [0.00784313725490196,0.505882352941176,0.541176470588235];[0.415686274509804,0.317647058823529,0.639215686274510];...
    [0.254901960784314,0.670588235294118,0.364705882352941]};
else 
    colors = ip.Results.colors; 
end 

if isempty(ip.Results.colorShadesSEM) 
    
colorShadesSEM = {[0.588235294117647,0.588235294117647,0.588235294117647]...
    [0.984313725490196,0.415686274509804,0.290196078431373]...
    [0.996078431372549,0.850980392156863,0.462745098039216]...
    [0.874509803921569,0.396078431372549,0.690196078431373]...
    [0.419607843137255,0.682352941176471,0.839215686274510]...
    [0.403921568627451,0.662745098039216,0.811764705882353]...
    [0.619607843137255,0.603921568627451,0.784313725490196]...
    [0.454901960784314,0.768627450980392,0.462745098039216]};
else 
    colorShadesSEM = ip.Results.colorShadesSEM; 
end 

colorShadesSD = {[0.850980392156863,0.850980392156863,0.850980392156863]...
[0.988235294117647,0.733333333333333,0.631372549019608]...
[1,1,0.800000000000000]...
[0.831372549019608,0.725490196078431,0.854901960784314]...
[0.776470588235294,0.858823529411765,0.937254901960784]...
[0.815686274509804,0.819607843137255,0.901960784313726]...
[0.854901960784314,0.854901960784314,0.921568627450980]...
[0.780392156862745,0.913725490196078,0.752941176470588]};
   
names = ip.Results.names; 
%% make NotBoxplot

valueMat = nanmean(dataPerCell,1); 
grpIDsVect = grpIDsVect'; 
FCondMean=figure;
hNB = notBoxPlot(valueMat,grpIDsVect);

nGroups = length(unique(grpIDsVect)); 
% Set the colors 
arrayfun(@(i) set(hNB(i).mu,'color',colors{i}),1:nGroups);
arrayfun(@(i) set(hNB(i).semPtch,'faceColor',colorShadesSEM{i}),1:nGroups);
arrayfun(@(i) set(hNB(i).sdPtch,'faceColor',colorShadesSD{i}(1,:)),1:nGroups);
arrayfun(@(i) set(hNB(i).data,'markerFaceColor', colors{i}),1:nGroups);
arrayfun(@(i) set(hNB(i).data,'markerEdgeColor','w'),1:nGroups);
%                 arrayfun(@(i) set(hNB(i).data,'markerSize',1),1:numel(names));
%                 arrayfun(@(i) uistack(hNB(i).sdPtch,'top'),1:numel(names));
%                 arrayfun(@(i) uistack(hNB(i).semPtch,'top'),1:numel(names));
%                 arrayfun(@(i) uistack(hNB(i).mu,'top'),1:numel(names));
if ~isempty(ip.Results.names)
    set(gca,'XTick',1:nGroups);
    set(gca,'XTickLabel',names,'FontSize',20);
    set(gca,'XTickLabelRotation',45);
end
h1 = get(gcf,'CurrentAxes');
yLim = h1.YLim; % use whatever they used
axis([0.5 2.5 yLim(1) yLim(2)]);
axis([0.5 nGroups + 0.5  yLim(1) yLim(2)]);
                
%% boxplot pooled 
if isempty(names)
    namesBoxplot = arrayfun(@(x) num2str(x),1:nGroups,'uniformoutput',0);
else
    namesBoxplot = names;
end 
FCondPooled=figure;
hBoxPool = boxplot(dataPerCell,grpIDsVect,'notch','on',...
                'symbol', '+','outliersize',1,'colors',vertcat(colors{1:nGroups}), ...
                'labelorientation','inline','Labels',namesBoxplot);
FPerCell=figure;
%% boxplot per cell 
hBoxPerCell=  boxplot(dataPerCell,'ColorGroup',grpIDsVect,'colors',vertcat(colors{1:nGroups}));
end

