function varargout = plotSkel(vertOrStruct,edgePaths,edgeLabels,mask,edgeCols)

%Um... add some documentation... sometime, like when you have some free
%time and would rather add documentation to your functions than do
%something fun or relaxing. Yeah, that would be a good time. Then.


if isstruct(vertOrStruct)
    vertices = vertOrStruct.vertices;
    edgePaths = vertOrStruct.edgePaths;
    edgeLabels = vertOrStruct.edgeLabels;
else
    vertices = vertOrStruct;
end


if nargout > 0
    figHan = fsFigure(.6);
end

if nargin > 3 && ~isempty(mask) 
    show3DMask(mask);
end

nEdges = numel(edgePaths);
if ~exist('edgeLabels','var') || isempty(edgeLabels)
    edgeLabels = ones(nEdges,1);
end

if nargin < 5 || isempty(edgeCols)
    edgeCols = jet(nEdges);
    singleCol = false;
elseif size(edgeCols,1) == 1
    singleCol = true;
    edgeCols = repmat(edgeCols,nEdges,1);
end

for j = 1:size(vertices,1)    
    if ~singleCol
        text(vertices(j,2),vertices(j,1),vertices(j,3),num2str(j),'Color','k')
    else
        text(vertices(j,2),vertices(j,1),vertices(j,3),num2str(j),'Color',edgeCols(1,:))
    end
    if j == 1
        hold on
    end
end


for j = 1:numel(edgePaths)
    if ~isempty(edgePaths{j})
        if edgeLabels(j) == 1
            plot3(edgePaths{j}(:,2),edgePaths{j}(:,1),edgePaths{j}(:,3),'color',edgeCols(j,:),'LineWidth',2);
        else
            plot3(edgePaths{j}(:,2),edgePaths{j}(:,1),edgePaths{j}(:,3),'--','color',edgeCols(j,:),'LineWidth',4);
        end
        text(mean(edgePaths{j}(:,2)),mean(edgePaths{j}(:,1)),mean(edgePaths{j}(:,3)),num2str(j),'color',edgeCols(j,:),'FontSize',12);
    end
end

view(3);
axis equal

if nargin < 4    
    set(gca,'Color',[.6 .6 .6]);
end

if nargout > 0
    varargout{1} = figHan;
end