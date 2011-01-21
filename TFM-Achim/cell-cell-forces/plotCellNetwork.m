function plotCellNetwork(network,xLimVal,yLimVal,scale)

edge  =network.edge;
node  =network.node;
maxMag=network.stats.maxMag;

if nargin < 4 || isempty(scale)
    scale=50/maxMag;
end

marker=['r','b','m','c','g','y','k'];
% plot the edges:
for j=1:length(edge)
    if ~isempty(edge{j})
        plot([edge{j}.strPt(1) edge{j}.endPt(1)],[edge{j}.strPt(2) edge{j}.endPt(2)],'-k')
        hold on
        if ~isempty(edge{j}.f1) && ~isempty(edge{j}.f2)
            % plot the first force, f1:
            quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),scale*edge{j}.f1(:,1),scale*edge{j}.f1(:,2),0,marker(mod(edge{j}.nodes(1),7)+1));
            % plot the second force, f2:
            quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),scale*edge{j}.f2(:,1),scale*edge{j}.f2(:,2),0,marker(mod(edge{j}.nodes(2),7)+1));
            plot(edge{j}.pos(:,1),edge{j}.pos(:,2),'.k','MarkerSize',10)
        else
            text(edge{j}.pos(:,1),edge{j}.pos(:,2),'???','VerticalAlignment','top','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))
        end
        % Use the results from the cluster Analysis if they exist
        if ~isempty(edge{j}.fc) % &&  isempty(edge{j}.f1) && isempty(edge{j}.f2)
            % plot the force obtained by the cluster analysis:
            quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),scale*edge{j}.fc(:,1),scale*edge{j}.fc(:,2),0,marker(mod(edge{j}.nodes(1),7)+1),'LineWidth',2);
            % plot the normal vector:
            quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),edge{j}.n_Vec(:,1),edge{j}.n_Vec(:,2),0,'k');
            plot(edge{j}.pos(:,1),edge{j}.pos(:,2),'.k','MarkerSize',10)
        end
        % plot the edge number
        text(edge{j}.pos(:,1),edge{j}.pos(:,2),[num2str(j),'^{',num2str(edge{j}.nodes),'}'],'VerticalAlignment','middle','HorizontalAlignment','center','color','r');% marker(mod(k,7)+1))
    end
end
% plot the residual forces of each cell:
for k=1:length(node)
    if ~isempty(node{k})
        quiver(node{k}.pos(:,1),node{k}.pos(:,2),node{k}.vec(:,1)*scale,node{k}.vec(:,2)*scale,0,'k') %marker(mod(k,7)+1))%,'MarkerSize',sizeCirc,'LineWidth',2)
        hold on
    end
end
% plot the number of each cell and its connectivity:
for k=1:length(node)
    if ~isempty(node{k})
        if isnan(node{k}.mag)
            sizeCirc=10^(-5);
        else
            sizeCirc=50*node{k}.mag/maxMag;
        end
        plot(node{k}.pos(:,1),node{k}.pos(:,2),['o',marker(mod(k,7)+1)],'MarkerFaceColor','w','MarkerSize',sizeCirc,'LineWidth',2)
        if isfield(node{k},'spec') && ~isempty(node{k}.spec) && node{k}.spec==1
            %text(node{k}.pos(:,1),node{k}.pos(:,2),['*',num2str(k),'^',num2str(node{k}.deg)],'VerticalAlignment','top','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))
            text(node{k}.pos(:,1),node{k}.pos(:,2),['*',num2str(k),'^',num2str(node{k}.deg),'_{',num2str(node{k}.edges),'}'],'VerticalAlignment','top','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))   
        else
            %text(node{k}.pos(:,1),node{k}.pos(:,2),[num2str(k),'^',num2str(node{k}.deg)],'VerticalAlignment','top','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))
            text(node{k}.pos(:,1),node{k}.pos(:,2),[num2str(k),'^',num2str(node{k}.deg),'_{',num2str(node{k}.edges),'}'],'VerticalAlignment','top','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))   
        end
        hold on
    end
end
% The scale bar for the stresses:
if nargin > 1 && ~isempty(xLimVal) && ~isempty(yLimVal)
    fxScaleBar_nN=100;
    fyScaleBar_nN=0;
    dPixX=130;
    dPixY=20;
    textSpace=20;
    quiver(xLimVal(2)-dPixX, yLimVal(2)-dPixY,scale*fxScaleBar_nN,scale*fyScaleBar_nN,0,'k','LineWidth',2,'MaxHeadSize',5)
    text(  xLimVal(2)-dPixX, yLimVal(2)-dPixY-textSpace,[num2str(fxScaleBar_nN),' nN'],'HorizontalAlignment','left','color', 'k','FontSize',16)
end
set(gca,'YDir','reverse')
% This is important for a proper scale bar!
axis equal
if nargin > 1 && ~isempty(xLimVal)
    xlim(xLimVal);
end
if nargin > 2 && ~isempty(yLimVal)
    ylim(yLimVal);
end
hold off

return;

% To plot and save the figures for one constrForceField:
padZeros=3;
for frame=1:31
    figure(frame)
    if ~isdir('cellNetwork')
        mkdir('cellNetwork')
    end
    filename = [pwd,filesep,'cellNetwork',filesep,'cellNetwork_',num2str(frame,['%0.',int2str(padZeros),'d'])];
    if isfield(constrForceField{frame},'network_tracked');
        plotCellNetwork(constrForceField{frame}.network_tracked,[350 1200],[300 750],0.5)
        saveas(gcf,[filename,'.tiff'],'tiffn');
%     saveas(gcf,[filename, '.eps'], 'psc2');
%     display(['Figure saved to: ',filename,'.tiffn+.eps'])
    end
end
    
