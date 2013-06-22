function [M]=calcFwdMap(x_vec_u, y_vec_u, forceMesh, E,span,meshPtsFwdSol,method)
if nargin<5
    baseSpan=1:forceMesh.numNodes;
    method = 'fast';
else
    baseSpan=span;
end

if nargin<6
    meshPtsFwdSol=[];
    method = 'fast';
end

ux=zeros(length(x_vec_u),2*forceMesh.numNodes);
uy=zeros(length(y_vec_u),2*forceMesh.numNodes);

for j=baseSpan
    display([num2str(j),' of: ',num2str(length(baseSpan))]);
    %limits for integration: integrate base function only over their
    %respective support:
    
    xmin=forceMesh.bounds(j).x(1);
    xmax=forceMesh.bounds(j).x(2);
    ymin=forceMesh.bounds(j).y(1);
    ymax=forceMesh.bounds(j).y(2);

    if strcmp(method,'fast')
        [ux(:,j),uy(:,j)]=fwdSolution(x_vec_u,y_vec_u,E,xmin,xmax,ymin,ymax,forceMesh.base(j).f_intp_x,forceMesh.base(j).f_intp_y,'fft',[],meshPtsFwdSol);
        [ux(:,j+forceMesh.numNodes),uy(:,j+forceMesh.numNodes)]=fwdSolution(x_vec_u,y_vec_u,E,xmin,xmax,ymin,ymax,forceMesh.base(j+forceMesh.numNodes).f_intp_x,forceMesh.base(j+forceMesh.numNodes).f_intp_y,'fft',[],meshPtsFwdSol);
    elseif strcmp(method,'conv_free')
        [ux(:,j),uy(:,j)]=fwdSolution(x_vec_u,y_vec_u,E,xmin,xmax,ymin,ymax,forceMesh.base(j).f_intp_x,forceMesh.base(j).f_intp_y,'conv_free',[],meshPtsFwdSol);
        [ux(:,j+forceMesh.numNodes),uy(:,j+forceMesh.numNodes)]=fwdSolution(x_vec_u,y_vec_u,E,xmin,xmax,ymin,ymax,forceMesh.base(j+forceMesh.numNodes).f_intp_x,forceMesh.base(j+forceMesh.numNodes).f_intp_y,'conv_free',[],meshPtsFwdSol);
    end
end

M=vertcat(ux,uy);

% plot an example to see if it works correctly
ind=10;
if forceMesh.numNodes>ind-1
    xmin=min(x_vec_u);
    ymin=min(y_vec_u);
    xmax=max(x_vec_u);
    ymax=max(y_vec_u);
    figure(11)
    quiver(x_vec_u,y_vec_u,ux(:,ind),uy(:,ind))
    hold on
    quiver(x_vec_u,y_vec_u,ux(:,ind+forceMesh.numNodes),uy(:,ind+forceMesh.numNodes))
    xlim([xmin xmax])
    ylim([ymin ymax])
    hold off
end