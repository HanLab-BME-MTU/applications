function plotComparisonInfo(in,conditions,outDir)
% Create box plots of the results of partitionComparison

plotter(in,conditions,'Fraction of tracks colocalizing',outDir,'fracInside',...
    false,'fracInside')
plotter(in,conditions,'Mean length of partitioned tracks',outDir,'meanLength',...
    false,'meanLengthInside','meanLengthOutside')
plotter(in,conditions,'Incidence of merge and splits per inside track',outDir,...
    'fracMergeSplitInside',false,'fracMergeInside','fracSplitInside')
plotter(in,conditions,'Incidence of merge and splits per outside track',outDir,...
    'fracMergeSplitOutside',false,'fracMergeOutside','fracSplitOutside')
plotter(in,conditions,'Asymmetry of tracks',outDir,'asymmetric',true,...
    'asymmetricInside','asymmetricOutside','notAsymmetricInside','notAsymmetricOutside')
plotter(in,conditions,'Diffusion classification using full dim.',outDir,...
    'classificationFullDim',true,'immobileFullInside','immobileFullOutside',...
    'confinedFullInside','confinedFullOutside','brownianFullInside',...
    'brownianFullOutside','directedFullInside','directedFullOutside')
plotter(in,conditions,'Diffusion classification using one dim.',outDir,...
    'classificationOneDim',true,'immobileOneInside','immobileOneOutside',...
    'confinedOneInside','confinedOneOutside','brownianOneInside',...
    'brownianOneOutside','directedOneInside','directedOneOutside')
plotter(in,conditions,...
    'Diffusion classification of pre-inside tracks using full dim.',outDir,...
    'classificationFullDimPrePost',true,'immobileFullPre','immobileFullPost',...
    'confinedFullPre','confinedFullPost','brownianFullPre','brownianFullPost',...
    'directedFullPre','directedFullPost')
plotter(in,conditions,...
    'Diffusion classification of post-inside tracks using full dim.',outDir,...
    'classificationOneDimPrePost',true,'immobileOnePre','immobileOnePost',...
    'confinedOnePre','confinedOnePost','brownianOnePre','brownianOnePost',...
    'directedOnePre','directedOnePost')
plotter(in,conditions,'Normalized diffusion coefficient with full dim.',outDir,...
    'normDiffCoefFull',false,'normDiffCoefFullInside','normDiffCoefFullOutside')
end

function plotter(info,conditions,plotTitle,outDir,outPath,showN,varargin)
% Make a box plot of the given data fields
% Inputs: showN displays number of data points (for fields that have that
% information stored in field "____N". varargin contains the fields to plot
nFields = numel(varargin);
data = zeros(0,1);
group = [];
n = {};
for iField = 1:nFields;
    for iCond = 1:numel(info)
        field = varargin{iField};
        newData = [info{iCond}.(field)];
        newData = newData(~isnan(newData));
        nData = numel(newData);
        data = [data;newData(:)];
        label = [field,'-',conditions{iCond}];
        group = [group;repmat({label},nData,1)];
        
        if showN && ~isempty(newData)
        	n = [n,{num2str(sum([info{iCond}.([field,'N'])]))}];
        end
    end
end
if isempty(data)
    % If no data, don't plot
    return
end
figure;
boxplot(data,group,'LabelOrientation','inline','Notch','on','Colors','mb')
title(plotTitle)
if showN
    for i = 1:numel(n)
        text(i-0.1,0,n{i},'FontSize',8)
    end
end
savefig([outDir,filesep,outPath])
close(gcf)
end




% function plotterBar(info,conditions,outDir,showN,varargin)
% data = [];
% error = [];
% n = {};
% for i = 1:numel(info)
%     for iField = 1:numel(varargin)
%         field = varargin{iField};
%         d = nanmean([info{i}.(field)]);
%         e = nanstd([info{i}.(field)]);
%         data(iField,i) = d;
%         error(iField,i) = e;
%         if showN
%         	n{iField,i} = num2str(sum([info{i}.([field,'N'])]));
%         end
%     end
% end
% figure
% h = bar(data);
% hold on
% 
% set(h,'BarWidth',1);    % The bars will now touch each other
% set(gca,'YGrid','on')
% set(gca,'GridLineStyle','-')
% 
% set(get(gca,'YLabel'),'String','U')
% 
% hold on;
% numgroups = numel(varargin); 
% numbars = numel(info); 
% groupwidth = min(0.8, numbars/(numbars+1.5));
% for i = 1:numbars
%       x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
%       errorbar(x, data(:,i), error(:,i), 'k.');
% end
% 
% 
% 
% % errorbar(1:numel(data),data,error,'.')
% labels = varargin;
% set(gca,'XTickLabel',labels);
% legend(gca,conditions);
% 
% 
% % if showN
% %     % label with n's
% %     xtick = get(gca,'xtick');
% %     nCon = numel(conditions);
% %     if mod(nCon,2) == 0
% %         offset = 1/nCon;
% %         xloc = [];
% %         center = (nCon+1)/2;
% %         for i = 1:nCon
% %             xloc = [xloc;xtick+2*offset*(i-center)];
% %         end
% %         xloc = xloc(:);
% %     else
% %         offset = (nCon-1)/4;
% %         xloc = [];
% %         center = (nCon+1)/2;
% %         for i = 1:nCon
% %             xloc = [xloc;xtick+offset*(i-center)];
% %         end
% %         xloc = xloc(:);
% %     end
% %     text(xloc,-0.5*ones(1,numel(n)),n{:})
% % end
% end errorbar(1:numel(data),data,error,'.')
% labels = varargin;
% set(gca,'XTickLabel',labels);
% legend(gca,conditions);
% 
% 
% % if showN
% %     % label with n's
% %     xtick = get(gca,'xtick');
% %     nCon = numel(conditions);
% %     if mod(nCon,2) == 0
% %         offset = 1/nCon;
% %         xloc = [];
% %         center = (nCon+1)/2;
% %         for i = 1:nCon
% %             xloc = [xloc;xtick+2*offset*(i-center)];
% %         end
% %         xloc = xloc(:);
% %     else
% %         offset = (nCon-1)/4;
% %         xloc = [];
% %         center = (nCon+1)/2;
% %         for i = 1:nCon
% %             xloc = [xloc;xtick+offset*(i-center)];
% %         end
% %         xloc = xloc(:);
% %     end
% %     text(xloc,-0.5*ones(1,numel(n)),n{:})
% % end
% end