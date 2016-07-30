function plotComparisonInfo(in,conditions,outDir)

plotter(in,conditions,[],false,'fracInside')
plotter(in,conditions,[],false,'meanLengthInside','meanLengthOutside')
% plotter(in,conditions,[],false,'fracMergeInside','fracMergeOutside','fracSplitInside','fracSplitOutside')
% plotter(in,conditions,[],true,'asymmetricInside','asymmetricOutside','notAsymmetricInside','notAsymmetricOutside')
% plotter(in,conditions,[],true,'immobileFullInside','immobileFullOutside','confinedFullInside','confinedFullOutside',...
%     'brownianFullInside','brownianFullOutside','directedFullInside','directedFullOutside')
% plotter(in,conditions,[],true,'immobileOneInside','immobileOneOutside','confinedOneInside','confinedOneOutside',...
%     'brownianOneInside','brownianOneOutside','directedOneInside','directedOneOutside')
% plotter(in,conditions,[],true,'immobileFullPre','immobileFullPost','confinedFullPre','confinedFullPost',...
%     'brownianFullPre','brownianFullPost','directedFullPre','directedFullPost')
% plotter(in,conditions,[],true,'immobileOnePre','immobileOnePost','confinedOnePre','confinedOnePost',...
%     'brownianOnePre','brownianOnePost','directedOnePre','directedOnePost')
end

function plotter(info,conditions,outDir,showN,varargin)
nFields = numel(varargin);
data = zeros(0,1);
group = [];
for iField = 1:nFields;
    for iCond = 1:numel(info)
        field = varargin{iField};
        newData = [info{iCond}.(field)];
        newData = newData(~isnan(newData));
        nData = numel(newData);
        data = [data;newData(:)];
        label = [field,'-',conditions{iCond}];
        group = [group;repmat({label},nData,1)];
    end
end
figure;
boxplot(data,group,'LabelOrientation','inline','Notch','marker')
end




function plotterBar(info,conditions,outDir,showN,varargin)
data = [];
error = [];
n = {};
for i = 1:numel(info)
    for iField = 1:numel(varargin)
        field = varargin{iField};
        d = nanmean([info{i}.(field)]);
        e = nanstd([info{i}.(field)]);
        data(iField,i) = d;
        error(iField,i) = e;
        if showN
        	n{iField,i} = num2str(sum([info{i}.([field,'N'])]));
        end
    end
end
figure
h = bar(data);
hold on

set(h,'BarWidth',1);    % The bars will now touch each other
set(gca,'YGrid','on')
set(gca,'GridLineStyle','-')

set(get(gca,'YLabel'),'String','U')

hold on;
numgroups = numel(varargin); 
numbars = numel(info); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
      x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
      errorbar(x, data(:,i), error(:,i), 'k.');
end



% errorbar(1:numel(data),data,error,'.')
labels = varargin;
set(gca,'XTickLabel',labels);
legend(gca,conditions);


% if showN
%     % label with n's
%     xtick = get(gca,'xtick');
%     nCon = numel(conditions);
%     if mod(nCon,2) == 0
%         offset = 1/nCon;
%         xloc = [];
%         center = (nCon+1)/2;
%         for i = 1:nCon
%             xloc = [xloc;xtick+2*offset*(i-center)];
%         end
%         xloc = xloc(:);
%     else
%         offset = (nCon-1)/4;
%         xloc = [];
%         center = (nCon+1)/2;
%         for i = 1:nCon
%             xloc = [xloc;xtick+offset*(i-center)];
%         end
%         xloc = xloc(:);
%     end
%     text(xloc,-0.5*ones(1,numel(n)),n{:})
% end
end