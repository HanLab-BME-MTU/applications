function [popRYGB,nPrctRYGB]=plusTipGroupDataQuadPlot(groupData,param1,thresh1,param2,thresh2,xLims,yLims)
% Makes quadrant scatterplot split by property percentiles or values, mapped onto image

% SYNOPSIS:
% [popRYGB,nPrctRYGB]=plusTipGroupDataQuadPlot(groupData,param1,thresh1,param2,thresh2,xLims,yLims)

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
%           thresh1/thresh2: values where data should be split for
%                            properties param1/param2, respectively.
%           xLims/yLims    : [lower upper] limits to x and y axes of
%                            scatter plot
%
% OUTPUT:
%           a scatterplot of the
%           parameters plotted against each other, split into four colors
%           at the boundaries created by the division marks
%
% 2009/08 - Kathryn Applegate Matlab R2008a
percentile1=[];
percentile2=[];
doPlot=1;

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
    msgbox('plusTipParamPlot: Input parameters must be from same group (i.e. growth)','Input Error','error')
    return
end

if nargin<6 || isempty(xLims)
    xLims=[];
end
if nargin<7 || isempty(yLims)
    yLims=[];
end


nGrps=length(groupData);
popRYGB=zeros(nGrps,4);
for iGrp=1:nGrps

    titleStr=groupData(iGrp,1).info.name;
    titleStr=strrep(titleStr,'_','-');

    % for param1 and param2, get data and label
    for iParam=1:2
        % assign data based on user input
        switch eval(['param' num2str(iParam)])

            case 'growthSpeed'
                tempData=groupData(iGrp,1).gs;
                tempLabel='Growth Speed (microns/min)';
            case 'fgapSpeed'
                tempData=groupData(iGrp,1).fs;
                tempLabel='Fgap Speed (microns/min)';
            case 'bgapSpeed'
                tempData=groupData(iGrp,1).bs;
                tempLabel='Bgap Speed (microns/min)';
            case 'growthLifetime'
                tempData=groupData(iGrp,1).gl;
                tempLabel='Growth Lifetime (sec)';
            case 'fgapLifetime'
                tempData=groupData(iGrp,1).fl;
                tempLabel='Fgap Lifetime (sec)';
            case 'bgapLifetime'
                tempData=groupData(iGrp,1).bl;
                tempLabel='Bgap Lifetime (sec)';
            case 'growthDisp'
                tempData=groupData(iGrp,1).gd;
                tempLabel='Growth Displacement (microns)';
            case 'fgapDisp'
                tempData=groupData(iGrp,1).fd;
                tempLabel='Fgap Displacement (microns)';
            case 'bgapDisp'
                tempData=groupData(iGrp,1).bd;
                tempLabel='Bgap Displacement (microns)';
            otherwise
                tempData='This parameter is not supported.';
        end

        % first iteration gives x-axis data, second gives y-axis data
        if iParam==1
            data1=tempData;
            if isempty(thresh1)
                thresh1=prctile(data1,percentile1);
            else
                percentile1=floor(100*length(find(data1<thresh1))/length(data1));
            end
            label1=[tempLabel ', P' num2str(percentile1) ' = ' sprintf('%3.2f',thresh1)];
        else
            data2=tempData;
            if isempty(thresh2)
                thresh2=prctile(data2,percentile2);
            else
                percentile2=floor(100*length(find(data2<thresh2))/length(data2));
            end
            label2=[tempLabel ', P' num2str(percentile2) ' = ' sprintf('%3.2f',thresh2)];
        end

    end

    colorMap=['b','g','y','r'];

    % make scatterplot showing data distribution between the two parameters
    if doPlot==1
        figure
    end
    for iColor=1:4
        switch iColor
            case 1 % blue
                idx=find(data1>thresh1 & data2>thresh2);
                popRYGB(iGrp,4)=length(idx);
            case 2 % green
                idx=find(data1<=thresh1 & data2>thresh2);
                popRYGB(iGrp,3)=length(idx);
            case 3 % yellow
                idx=find(data1>thresh1 & data2<=thresh2);
                popRYGB(iGrp,2)=length(idx);
            case 4 % red
                idx=find(data1<=thresh1 & data2<=thresh2);
                popRYGB(iGrp,1)=length(idx);
        end
        if doPlot==1
            %         if strcmp(param1,'growthSpeed') && strcmp(param2,'growthLifetime') && ~isempty(xLims) && ~isempty(yLims)
            %             [x,y]=meshgrid(1:xLims(2),1:yLims(2));
            %             disp=x.*(y./60); % displacement in microns
            %             contourf(disp)
            %             colormap gray
            %             hold on
            %         end
            scatter(data1(idx),data2(idx),[],colorMap(iColor),'.');
            hold on
        end
    end
    if doPlot==1
        if ~isempty(xLims)
            xlim(xLims)
        end
        if ~isempty(yLims)
            ylim(yLims)
        end

        xlabel(label1)
        ylabel(label2)

        % any track which ends before the frame range begins or begins after the frame range ends will be excluded.
        title([titleStr ', N = ' num2str(length(data1)) ' tracks']);
    end
end
[nPrctRYGB]=plusTipQuadColorbar(popRYGB);

