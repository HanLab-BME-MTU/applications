E=100;  %Young's modulus
L=0;
meshPtsFwdSol=2^10;

numPoints_u=10       %(muss gerade Anzahl sein)
numPoints_f=10       %(muss gerade Anzahl sein)
numPoints_out=20

xmin =-5;
xmax = 5;
ymin =-5;
ymax = 5;

[x_mat_u y_mat_u]=meshgrid(linspace(xmin,xmax,numPoints_u) , linspace(ymin,ymax,numPoints_u));
x_vec_u=reshape(x_mat_u,[],1);
y_vec_u=reshape(y_mat_u,[],1);

[x_mat_f y_mat_f]=meshgrid(linspace(xmin,xmax,numPoints_f) , linspace(ymin,ymax,numPoints_f));
x_vec_f=reshape(x_mat_f,[],1);
y_vec_f=reshape(y_mat_f,[],1);

[ux uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,@(x,y) assumedForce(1,x,y),@(x,y) assumedForce(2,x,y),'fft',[],meshPtsFwdSol)
ux_vec=reshape(ux,[],1);
uy_vec=reshape(uy,[],1);
u=vertcat(ux_vec,uy_vec);

figure(1)
quiver(x_mat_u,y_mat_u,assumedForce(1,x_mat_u,y_mat_u),assumedForce(2,x_mat_u,y_mat_u));

figure(2)
quiver(x_mat_u,y_mat_u,ux,uy);


%***************** here starts the BEM-reconstruction *********************
['expected computation time:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)),'s these are:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)/3600),'h']    

forceMesh=createMeshAndBasis(x_vec_f,y_vec_f);

nTasks=ceil(numPoints_u^2*numPoints_f^2*27.6*10^(-3)/(1.5*60));
spanTask=floor(forceMesh.numNodes/nTasks);

k=1;
toDoList=1:forceMesh.numNodes;
while ~isempty(toDoList)
    if length(toDoList)>spanTask
        imputFile.task(k).span=toDoList(1:spanTask);
        toDoList(1:spanTask)=[];
    else
        imputFile.task(k).span=toDoList(1:end);
        toDoList=[];
    end
    k=k+1;
end

imputFile.E=E;
imputFile.L=L;
imputFile.x=x_vec_u;
imputFile.y=y_vec_u;
imputFile.u=u;
imputFile.numTasks=length(imputFile.task);
imputFile.forceMesh=forceMesh;


save([pwd,'/files/imputFile.mat'],'imputFile');


