function [displField, forceField, out]=createDisplField(sdT,inputFileList,target_dir,filter,yModu_Pa,pRatio,pixSize_mu,regParam,method,meshPtsFwdSol,xrange,yrange,doRotReg)
saveAllBEMpar=1;

if nargin < 1 || isempty(sdT)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select ResidualT.mat to be used');
       %the stage drift Transformation:
       fileStruct=load([pathname filesep filename]);
       sdT=fileStruct.T;
elseif ischar(sdT)
       fileStruct=load(sdT);
       sdT=fileStruct.T;
end

%read in stack of flow files:
if nargin < 2 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select the first flow file to be analyzed');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   inputFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for i = 1:numel(inputFileList)
        isValid = isValid && exist(inputFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

%get the target directory:
if nargin < 3 || isempty(target_dir)
    targetDir = uigetdir('','Select target directory for displField.mat');
elseif ischar('target_dir')
    targetDir=target_dir;
end

if nargin < 4 || isempty(filter)
    filter=[];
elseif length(filter)==1
    numStd=filter(1);
    boxSizeLocFac=[];
    boxSizeGlbFac=[];
    maxItr=1;
elseif length(filter)==2
    numStd=filter(1);
    boxSizeLocFac=filter(2);
    boxSizeGlbFac=[];
    maxItr=1;
elseif length(filter)==3
    numStd=filter(1);
    boxSizeLocFac=filter(2);
    boxSizeGlbFac=filter(3);
    maxItr=1;
elseif length(filter)==4
    numStd=filter(1);
    boxSizeLocFac=filter(2);
    boxSizeGlbFac=filter(3);
    maxItr=filter(4);
end

if nargin < 5 || isempty(yModu_Pa)
    yModu_Pa=20000;
end

if nargin < 6 || isempty(pRatio)
    pRatio=0.5;
end

if nargin < 7 || isempty(pixSize_mu)
    pixSize_mu=0.163;
end

if nargin < 8 || isempty(regParam)
    regParam=10^(-7);
end

if nargin < 9 || isempty(method)
    method='FTTC';
end

if nargin < 10 || isempty(meshPtsFwdSol)
    meshPtsFwdSol=2^11;
end

if nargin < 13  || isempty(doRotReg)
    doRotReg=0;
end


% if the field of view is too large one has the option to choose a
% x-y-range where the force field should be calculated.
if nargin <11 || isempty(xrange) || isempty(yrange)
    xrange=[];
    yrange=[];
end

if isempty(filter)
    doPlot=0;
else
    doPlot=1;
end

display('Used parameters are:');
display(['filter: ',num2str(filter)]);
display(['yModu_Pa: ',num2str(yModu_Pa)]);
display(['pRatio: ',num2str(pRatio)]);
display(['pixSize_mu: ',num2str(pixSize_mu)]);
display(['regParam: ',num2str(regParam)]);
display(['meshPtsFwdSol: ',num2str(meshPtsFwdSol)]);



n = numel(inputFileList);

out.vec=[];
out.pos=[];
out.num=[];
out.filter=filter;

for i=1:n
    text='Status';
    progressText(i/n,text);
    fileStruct=load(inputFileList{i});
    flow=fileStruct.flow;
    flow(isnan(flow(:,3)),:)=[];
    flow(isnan(flow(:,4)),:)=[];
    
    %If available, perform stage drift correction:
    if ~isempty(sdT)
        %if length(sdT)==n maybe +-1
        Tx=sdT(i,1);
        Ty=sdT(i,2);    
        %else            
            %error('Size of the stage drift correction transform doesnt match the number of flow files');
        %end
    else
        Tx=0;
        Ty=0;
    end
    
    % If a crop region is chosen, cut out the relevant displacementfield:
    if ~isempty(xrange)
        flow(flow(:,2)<xrange(1),:)=[];
        flow(flow(:,2)>xrange(2),:)=[];
    end
    
    if ~isempty(yrange)
        flow(flow(:,1)<yrange(1),:)=[];
        flow(flow(:,1)>yrange(2),:)=[];
    end
    
    
    %the +sign in front of Ty is correct and also the +Tx.
    displField(i).vec=[flow(:,4)-flow(:,2)+Ty flow(:,3)-flow(:,1)+Tx];
        
    %The coordinates in the reference frame must not be transformed:
    displField(i).pos=[flow(:,2) flow(:,1)];
    
    displField(i).par.filter=filter;
   
 
    if ~isempty(filter)
        pos_fltr=displField(i).pos;
        vec_fltr=displField(i).vec;
        pos_out_all=[];
        vec_out_all=[];
        id_out_all =[];
        for iItr=1:maxItr
            [pos_fltr,vec_fltr,pos_out,vec_out,~,id_out]=filterVectorOutliers(pos_fltr,vec_fltr,numStd,boxSizeLocFac,boxSizeGlbFac,[]);
            pos_out_all=vertcat(pos_out_all,pos_out);
            vec_out_all=vertcat(vec_out_all,vec_out);
            % the actual id's are wrong when maxIt>1;
            id_out_all=vertcat(id_out_all,id_out);            
        end
        displField(i).pos = pos_fltr;
        displField(i).vec = vec_fltr;
        out(i).pos = pos_out_all;
        out(i).vec = vec_out_all;
        out(i).num = length(id_out_all);
            
        % dummy=input('Press enter to proceed: ');
%         figure()
%         quiver(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1),displField(i).vec(:,2),0)
%         hold on;
%         quiver(out(i).pos(:,1),out(i).pos(:,2),out(i).vec(:,1),out(i).vec(:,2),0,'r')
%         hold off;

    end
    displField(i).par.prep4fastBEM=0;
end
display('Number of points filtered out: '); 
display([num2str((1:n)'),repmat(': ',n,1),num2str([out(:).num]')]);

if doRotReg
   displField=perfRotReg(displField,1);
end

if strcmp(method,'FastBEM')
   displField=prepDisplForBEM(displField,'linear');
end

save([targetDir,filesep,'displField.mat'], 'displField');
save([targetDir,filesep,'out.mat'], 'out');

for i=1:length(displField)
    maxX(i)=max(displField(i).pos(:,1));
    maxY(i)=max(displField(i).pos(:,2));
end
theXlim=max(maxX);
theYlim=max(maxY);

padZeros=floor(log10(length(displField)))+1;

targetDirDispl=[targetDir,filesep ,'displField', filesep];
if ~isdir(targetDirDispl)
    mkdir(targetDirDispl)
end

if doPlot
    for i=1:n
        figure(1)
        gridSize=ceil(sqrt((theXlim)*(theYlim)/length(displField(i).vec)));
        maxForcePlot=0.3/gridSize*max(sqrt(displField(i).vec(:,1).^2+displField(i).vec(:,2).^2));

        quiver(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1)/maxForcePlot,displField(i).vec(:,2)/maxForcePlot,0,'k')
        if ~isempty(filter)
            hold on
            quiver(out(i).pos(:,1),out(i).pos(:,2),out(i).vec(:,1)/maxForcePlot,out(i).vec(:,2)/maxForcePlot,0,'r')
            hold off
        end
        xlim([1 theXlim])
        ylim([1 theYlim])
        set(gca,'DataAspectRatio', [1,1,50],'YDir','reverse')%,'XTick',[],'YTick',[])
        title(['Displacement field frame no: ',num2str(i),'. Outliers are red!'])
        if ~isempty(filter)
            if out(i).num==0
                text(theXlim/2,theYlim/2,'No Outliers');
            end
        end       
        saveas(gcf,[targetDirDispl, 'displField',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
    end
end

%*****************************************
% here begins the force reconstruction   *
%*****************************************

% For Benedikt's software to work, the displacement field has to be
% interpolated on a rectangular grid, with an even number of grid points
% along each edge. Furthermore, one has to make sure that the noisy data 
% has not to be extrapolated. This may happen along the edges. To prevent
% this, extract the corner of the displacement grid, calculate how often 
% (even number) the optimal gridspacing fits into each dimension, then 
% place the regular grid centered to the orignal bounds. Thereby make sure 
% that the edges have been eroded to a certain extend. This is performed by
% the following function.

[reg_grid,~,~,gridSpacing]=createRegGridFromDisplField(displField);


targetDirForce=[targetDir,filesep ,'forceField', filesep];
if ~isdir(targetDirForce)
    mkdir(targetDirForce)
end

for i=1:length(displField)
    [grid_mat,iu_mat, i_max,j_max] = interp_vec2grid(displField(i).pos, displField(i).vec,[],reg_grid);
    
    if strcmp(method,'FastBEM')
        % If grid_mat=[], then an optimal hexagonal force mesh is created
        % given the bead locations defined in displField:
        tic;
        if i==1 || displField(i).par.prep4fastBEM==0;
            [pos_f, force, forceMesh, M, pos_u, u, sol_coef, sol_mats]=reg_FastBEM_TFM(grid_mat, displField, i, yModu_Pa, pRatio, regParam, meshPtsFwdSol);
            display('The total time for calculating the FastBEM solution: ')
        elseif i>1 && displField(i).par.prep4fastBEM==1
            % since the displ field has been prepared such
            % that the measurements in different frames are ordered in the
            % same way, we don't need the position information any
            % more. The displ. measurements are enough.
            display('5.) Re-evaluate the solution:... ')
            % pull the new u-vector:
            u=vertcat(displField(i).vec(:,1),displField(i).vec(:,2));
            % recalculate the solution for the new displacement vec:
            [pos_f,force,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,u,forceMesh,regParam,[],[]);
            display(['Done: solution for frame: ',num2str(i)]);
        end
        toc;
        
        % The following values should/could be stored for the BEM-method.
        % In most cases, except the sol_coef this has to be stored only
        % once for all frames!
        if saveAllBEMpar==1
            forceField(i).par.forceMesh     = forceMesh;
            forceField(i).par.sol_coef      = sol_coef;
            forceField(i).par.M             = M; % This should not be saved every time! Although necessary to calculate the L-curve!
            forceField(i).par.sol_mats      = sol_mats;
            forceField(i).par.pos           = pos_u;
            forceField(i).par.u             = u;  
            forceField(i).par.meshPtsFwdSol = meshPtsFwdSol;   
        end
    else
        [pos_f,~,force,~,~,~] = reg_fourier_TFM(grid_mat, iu_mat, yModu_Pa, pRatio, pixSize_mu, gridSpacing, i_max, j_max, regParam);
    end   
    
    if doPlot
        figure(100)
        quiver(pos_f(:,1),pos_f(:,2),force(:,1),force(:,2))
        xlim([1 theXlim])
        ylim([1 theYlim])
    %    set(gca,'YDir','reverse')%,'XTick',[],'YTick',[])
        title(['Force field frame no: ',num2str(i)])
        saveas(gcf,[targetDirForce,'forceField',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
        set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse')%,'XTick',[],'YTick',[])
    end
    
    
    % Fill in the values to be stored:
    forceField(i).pos=pos_f;
    forceField(i).vec=force;
    forceField(i).posShifted=[]; % this will be calculated when needed
    forceField(i).vecReIntp =[]; % this will be calculated when needed
    forceField(i).par.yModu_Pa   =yModu_Pa;
    forceField(i).par.pRatio     =pRatio;
    forceField(i).par.pixSize_mu =pixSize_mu;
    forceField(i).par.regParam   =regParam;
    forceField(i).par.gridSpacing=gridSpacing;
    forceField(i).par.method     =method;
    forceField(i).doc='Correct pairings: (pos,vec):raw, (posShifted,vec):shifted according to displ, (pos,vecReIntp):reinterp to original grid';
    
    clear grid_mat;
    clear iu;
    clear iu_mat;
    
    %save([targetDir,filesep, 'forceField.mat'], 'forceField', '-v7.3');
end
display('Saving the forceField. In case of "FastBEM" this might take some time since the forward map is huge:...')
save([targetDir,filesep, 'forceField.mat'], 'forceField', '-v7.3');