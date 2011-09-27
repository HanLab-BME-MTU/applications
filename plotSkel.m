function varargout = plotSkel(vertices,edgePaths,edgeLabels,mask)


if nargout > 0
    figHan = fsFigure(.6);
end

if nargin > 3 && ~isempty(mask) 
    show3DMask(mask);
end

for j = 1:size(vertices,1)    
    text(vertices(j,2),vertices(j,1),vertices(j,3),num2str(j),'Color','k')
    if j == 1
        hold on
    end
end

nEdges = numel(edgePaths);
edgeCols = jet(nEdges);

if nargin < 3 || isempty(edgeLabels)
    edgeLabels = ones(nEdges,1);
end

for j = 1:numel(edgePaths)
    if ~isempty(edgePaths{j})
        if edgeLabels(j) == 1
            plot3(edgePaths{j}(:,2),edgePaths{j}(:,1),edgePaths{j}(:,3),'color',edgeCols(j,:),'LineWidth',2);
        else
            plot3(edgePaths{j}(:,2),edgePaths{j}(:,1),edgePaths{j}(:,3),'--','color',edgeCols(j,:),'LineWidth',4);
        end
        text(mean(edgePaths{j}(:,2)),mean(edgePaths{j}(:,1)),mean(edgePaths{j}(:,3)),num2str(j),'color',edgeCols(j,:));
    end
end

view(3);
axis equal

if nargout > 0
    varargout{1} = figHan;
end