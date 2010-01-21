function plusTipMakeBoxplots(data,setName,saveDir)
% plusTipMakeBoxplots saves speed, lifetime, and displacement boxplots
% for growth, fgap, and bgap populations from different movie groups
%
% SYNOPSIS: plusTipMakeBoxplots(data,setName,saveDir)
%
% see plusTipPoolGroupData



if ~isdir(saveDir)
    mkdir(saveDir)
end

% make within-group boxplots
gIdxInd=cellfun(@(x) find(x(:,5)==1),data,'UniformOutput',0);
fIdxInd=cellfun(@(x) find(x(:,5)==2),data,'UniformOutput',0);
bIdxInd=cellfun(@(x) find(x(:,5)==3),data,'UniformOutput',0);
for iParam=1:3
    switch iParam
        case 1
            if isempty(cell2mat(gIdxInd))
                continue
            end
            ms=cellfun(@(x,y) x(y,4),data,gIdxInd,'UniformOutput',0)';
            ml=cellfun(@(x,y) x(y,6),data,gIdxInd,'UniformOutput',0)';
            md=cellfun(@(x,y) x(y,7),data,gIdxInd,'UniformOutput',0)';
            titleStr='growth';
        case 2
            if isempty(cell2mat(fIdxInd))
                continue
            end
            ms=cellfun(@(x,y) x(y,4),data,fIdxInd,'UniformOutput',0)';
            ml=cellfun(@(x,y) x(y,6),data,fIdxInd,'UniformOutput',0)';
            md=cellfun(@(x,y) x(y,7),data,fIdxInd,'UniformOutput',0)';
            titleStr='fgap';
        case 3
            if isempty(cell2mat(bIdxInd))
                continue
            end
            ms=cellfun(@(x,y) x(y,4),data,bIdxInd,'UniformOutput',0)';
            ml=cellfun(@(x,y) x(y,6),data,bIdxInd,'UniformOutput',0)';
            md=cellfun(@(x,y) x(y,7),data,bIdxInd,'UniformOutput',0)';
            titleStr='bgap';
    end
    maxSize=cellfun(@(x) length(x),ms);
    allMS=nan(max(maxSize,[],2),length(ms));
    for i=1:length(ms)
        allMS(1:maxSize(i),i) = ms{1,i};
    end
    allML=nan(max(maxSize,[],2),length(ml));
    for i=1:length(ml)
        allML(1:maxSize(i),i) = ml{1,i};
    end
    allMD=nan(max(maxSize,[],2),length(md));
    for i=1:length(md)
        allMD(1:maxSize(i),i) = md{1,i};
    end
    % for each value, put group name into matrix
    grpVar=repmat(setName',[max(maxSize'),1]);

    figure % boxplot of speeds
    boxplot(allMS(:),grpVar(:),'notch','on','orientation','horizontal');
    title([titleStr ' speeds'])
    set(gca,'YDir','reverse')
    xlabel('speed (microns/min)')
    saveas(gcf,[saveDir filesep 'boxplot_speed_' titleStr '.fig'])
    saveas(gcf,[saveDir filesep 'boxplot_speed_' titleStr '.tif'])
    close(gcf)

    figure % boxplot of lifetimes
    boxplot(allML(:),grpVar(:),'notch','on','orientation','horizontal');
    title([titleStr ' lifetimes'])
    set(gca,'YDir','reverse')
    xlabel('lifetimes (sec)')
    saveas(gcf,[saveDir filesep 'boxplot_lifetime_' titleStr '.fig'])
    saveas(gcf,[saveDir filesep 'boxplot_lifetime_' titleStr '.tif'])
    close(gcf)

    figure % boxplot of displacements
    boxplot(allMD(:),grpVar(:),'notch','on','orientation','horizontal');
    title([titleStr ' displacements'])
    set(gca,'YDir','reverse')
    xlabel('displacement (microns)')
    saveas(gcf,[saveDir filesep 'boxplot_displacement_' titleStr '.fig'])
    saveas(gcf,[saveDir filesep 'boxplot_displacement_' titleStr '.tif'])
    close(gcf)
end