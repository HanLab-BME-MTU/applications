function plotCellNetwork(network,xLimVal,yLimVal,scale,noErrs,noTicks,noFn)
% for frame=1:length(network)
%     if ~isempty(network{frame})
%         display(['Plot frame: ',num2str(frame)]);
%         plotCellNetwork(network{frame})
%         pause(1)
%     end
% end

if nargin<5 || isempty(noErrs) || noErrs
    if network.stats.errs>0
        display('Found at least one node with unresolved force field! To plot this use noErrs=0')
        return
    end
end

edge  =network.edge;
node  =network.node;
maxMag=network.stats.maxMag;

if nargin < 4 || isempty(scale)
    scale=50/maxMag;
end

if nargin < 6 || isempty(noTicks)
    noTicks = 0;
end

if nargin < 7 || isempty(noFn)
    noFn = 0;
end


marker=['r','b','m','c','g','y','k'];
% plot the edges:
for j=1:length(edge)
    if ~isempty(edge{j})
        plot([edge{j}.strPt(1) edge{j}.endPt(1)],[edge{j}.strPt(2) edge{j}.endPt(2)],'-k')
        hold on
        if ~noFn && ~isempty(edge{j}.f1) && ~isempty(edge{j}.f2)
            % plot the first force, f1:
            quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),scale*edge{j}.f1(:,1),scale*edge{j}.f1(:,2),0,marker(mod(edge{j}.nodes(1),7)+1));
            % plot the second force, f2:
            quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),scale*edge{j}.f2(:,1),scale*edge{j}.f2(:,2),0,marker(mod(edge{j}.nodes(2),7)+1));
            plot(edge{j}.pos(:,1),edge{j}.pos(:,2),'.k','MarkerSize',10)
        elseif ~noFn
            text(edge{j}.pos(:,1),edge{j}.pos(:,2),'???','VerticalAlignment','top','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))
        end
        % Use the results from the cluster Analysis if they exist
        if ~isnan(sum(edge{j}.fc)) % &&  isempty(edge{j}.f1) && isempty(edge{j}.f2)
            % plot the force obtained by the cluster analysis:
            quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),scale*edge{j}.fc1(:,1),scale*edge{j}.fc1(:,2),0,marker(mod(edge{j}.nodes(1),7)+1),'LineWidth',1);
            % quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),scale*edge{j}.fc(:,1) ,scale*edge{j}.fc(:,2) ,0,marker(mod(edge{j}.nodes(1),7)+1),'LineWidth',2);
            % plot the normal vector:
            quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),scale*edge{j}.fc2(:,1),scale*edge{j}.fc2(:,2),0,marker(mod(edge{j}.nodes(2),7)+1),'LineWidth',1);
            quiver(edge{j}.pos(:,1),edge{j}.pos(:,2),edge{j}.n_Vec(:,1),edge{j}.n_Vec(:,2),0,'k');
            plot(edge{j}.pos(:,1),edge{j}.pos(:,2),'.k','MarkerSize',10)
        end
        % find the pair degree:
        conNodes=edge{j}.nodes;
        pairDeg=sort([node{conNodes(1)}.deg,node{conNodes(2)}.deg]);        
        % plot the edge number with id of connected cells:
        % text(edge{j}.pos(:,1),edge{j}.pos(:,2),[num2str(j),'^{',num2str(edge{j}.nodes),'}'],'VerticalAlignment','middle','HorizontalAlignment','center','color','r');% marker(mod(k,7)+1))
        % plot the edge number with the pair-degree:
        text(edge{j}.pos(:,1),edge{j}.pos(:,2),[num2str(j),'^{',num2str(pairDeg),'}'],'VerticalAlignment','middle','HorizontalAlignment','center','color','r');% marker(mod(k,7)+1))
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
        % circle size equals the magnitude of the residual force:
        % if isnan(node{k}.mag)
        %    sizeCirc=10^(-5);
        % else
        %    sizeCirc=50*node{k}.mag/maxMag;
        % end
        
        % circle size equals the sum of all interface forces (magnitude):
        sumIntF=sumIntForces(k,node,edge);
        if isnan(sumIntF)
            sizeCirc=10^(-5);
        else
            sizeCirc=13*sumIntF/maxMag;
        end
        if isfield(node{k},'spec') && ~isempty(node{k}.spec) && node{k}.spec==1
            plot(node{k}.pos(:,1),node{k}.pos(:,2),['o',marker(mod(k,7)+1)],'MarkerFaceColor','r','MarkerSize',sizeCirc,'LineWidth',2)
            text(node{k}.pos(:,1),node{k}.pos(:,2),['*',num2str(k),'^',num2str(node{k}.deg)],'VerticalAlignment','middle','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))
            %text(node{k}.pos(:,1),node{k}.pos(:,2),['*',num2str(k),'^',num2str(node{k}.deg),'_{',num2str(node{k}.edges),'}'],'VerticalAlignment','top','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))   
        else
            plot(node{k}.pos(:,1),node{k}.pos(:,2),['o',marker(mod(k,7)+1)],'MarkerFaceColor','w','MarkerSize',sizeCirc,'LineWidth',2)
            text(node{k}.pos(:,1),node{k}.pos(:,2),[num2str(k),'^',num2str(node{k}.deg)],'VerticalAlignment','middle','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))
            %text(node{k}.pos(:,1),node{k}.pos(:,2),[num2str(k),'^',num2str(node{k}.deg),'_{',num2str(node{k}.edges),'}'],'VerticalAlignment','top','HorizontalAlignment','center','color','k');% marker(mod(k,7)+1))   
        end
        hold on
    end
end
% The scale bar for the stresses:
if nargin > 1 && ~isempty(xLimVal) && ~isempty(yLimVal)
    fxScaleBar_nN=300;
    fyScaleBar_nN=0;
    dPixX=130;
    dPixY=20;
    textSpace=20;
    % lower right corner:
    % quiver(xLimVal(2)-dPixX, yLimVal(2)-dPixY,scale*fxScaleBar_nN,scale*fyScaleBar_nN,0,'k','LineWidth',2,'MaxHeadSize',5)
    % text(  xLimVal(2)-dPixX, yLimVal(2)-dPixY-textSpace,[num2str(fxScaleBar_nN),' nN'],'HorizontalAlignment','left','color', 'k','FontSize',16)
    % upper right corner:
    quiver(xLimVal(2)-1.25*dPixX, yLimVal(1)+4*dPixY,scale*fxScaleBar_nN,scale*fyScaleBar_nN,0,'k','LineWidth',2,'MaxHeadSize',5)
    text(  xLimVal(2)-1.25*dPixX, yLimVal(1)+4*dPixY-textSpace,[num2str(fxScaleBar_nN),' nN'],'HorizontalAlignment','left','color', 'k','FontSize',16)
end
if noTicks
    set(gca,'YDir','reverse','XTick',[],'YTick',[]);
else
    set(gca,'YDir','reverse');
end
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

function sumF=sumIntForces(nodeIdFunc,nodeFunc,edgeFunc)
%find all edges connected to that node:
edgesFunc=nodeFunc{nodeIdFunc}.edges;

%get the interface forces that belong to that node:
sumF=0;
for edgeIdFunc=edgesFunc
    if ~isempty(edgeFunc{edgeIdFunc}.f1) && ~isempty(edgeFunc{edgeIdFunc}.f2) && ~isnan(edgeFunc{edgeIdFunc}.f1(1)) && ~isnan(edgeFunc{edgeIdFunc}.f2(1))
        % then take the network force:
        fEdge=1/2*(edgeFunc{edgeIdFunc}.f1-edgeFunc{edgeIdFunc}.f2);
    else
        fEdge=edgeFunc{edgeIdFunc}.fc;
    end
    sumF=sumF+sqrt(sum(fEdge.^2));
end
    
