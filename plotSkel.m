function figHan = plotSkel(vertices,edges,edgePaths)

figHan = fsFigure(.75);
hold on

for j = 1:size(vertices,1)
    %plot3(vertices(j,2),vertices(j,1),vertices(j,3),'ok','MarkerSize',15)
    text(vertices(j,2),vertices(j,1),vertices(j,3),num2str(j),'Color','k')
end

nEdges = numel(edgePaths);
edgeCols = jet(nEdges);

for j = 1:numel(edgePaths)
    plot3(edgePaths{j}(:,2),edgePaths{j}(:,1),edgePaths{j}(:,3),'color',edgeCols(j,:));
    text(mean(edgePaths{j}(:,2)),mean(edgePaths{j}(:,1)),mean(edgePaths{j}(:,3)),num2str(j),'color',edgeCols(j,:));
end

view(3);
axis equal