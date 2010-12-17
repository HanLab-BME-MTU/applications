function [constrForceField]=cutOutForceFieldManyCells(forceField,imageFileList,path_cellCellForces,doIntp)

if nargin < 1 || isempty(forceField) 
   [filename_forceField, pathname_forceField] = uigetfile({'*.mat';'*.*'}, ...
       'Select forceField.mat to be used');
       %the stage drift Transformation:
       fileStruct=load([pathname_forceField filesep filename_forceField]);
       forceField=fileStruct.forceField;
       
       %target_dir=pathname_forceField;
end

%read in Stack of images:
if nargin < 2 || isempty(imageFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Ecad-Image');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   imageFileList = getFileStackNames([pathname filesep filename]);
elseif isdir(imageFileList)
    imageFileList=getFileListFromFolder(imageFileList);
else
    isValid = 1;
    for i = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end


if nargin < 4 || isempty(doIntp)
    doIntp=1;
    display('The interpolation will be performed at the very end of the script!')
    display('If the program crashes you simply have to execute: updateConstrForceField which performs the interpolation!')
end
%method='noIntp';

finished=false;
constrForceField=[];
toDoList=1:length(forceField);
numCells=NaN;
numCellsOld=-1;

while finished==false    
    try
        %file_path=[target_dir,filesep,'cellCellForces.mat'];
        load(path_cellCellForces);
        display(['Done files are: 1...',num2str(length(constrForceField)),' of in total: 1...',num2str(length(forceField))])
        reply=input('Which fields do you want to re-analyze, Press Enter for All, or name frames as [n1 n2 ...]: ');
        if isempty(reply)
            toDoList=1:length(forceField);
            for i=toDoList
                if i<=length(constrForceField)
                    constrForceField{i}=[];
                end
            end
        else
            toDoList=reply;
            for i=toDoList
                if i<=length(constrForceField)
                    constrForceField{i}=[];
                end
            end
        end
    catch exception
        if strcmp(exception.identifier, 'MATLAB:load:couldNotReadFile')
            display('No "cellCellForces.mat" file found, run through analysis completely.')
        else
            display(['Unknown error: ',exception.identifier])
        end
    end

    dilationR=35;
    sigmaGauss=[];
    closureRadius=[];
    % The following values are optional:
    holes=1;
    thrHoleLen=50;
    
    
    numFrames=length(toDoList);
    saveInterval=min(5,ceil(0.1*numFrames));

    %toDoList=30;
    for i=toDoList
        display('Plotting constraint forceField might take some time...')
        currentImage = double(imread(imageFileList{i}));
        roiOK='n';
        while strcmp(roiOK,'n') || strcmp(roiOK,'no')
            cellEdgeResults=cellPerim(imageFileList{i},dilationR,sigmaGauss,closureRadius,holes,thrHoleLen);
            segmRes{i}=cellEdgeResults{1};
            
            % store the segmentation results:
            constrForceField{i}.segmRes=segmRes{i};            
            constrForceField{i}.par=forceField(i).par;
            ImgSize=size(constrForceField{i}.segmRes.maskDilated);
            
            fpos=forceField(i).pos;
            fvec=forceField(i).vec;
            if sum(sum(abs(fpos-round(fpos))))==0
                [inpos,invec]=findVectorFieldInMask(fpos,fvec,constrForceField{i}.segmRes.maskDilated);
                
                constrForceField{i}.roi.pos=inpos;
                constrForceField{i}.roi.vec=invec;
            else
                display(['Use inpolygon for frame: ',num2str(i)]);
                checkVector = inpolygon(fpos(:,1),fpos(:,2),constrForceField{i}.segmRes.curveDilated(:,1),constrForceField{i}.segmRes.curveDilated(:,2));
                
                constrForceField{i}.roi.pos=fpos(checkVector,:);
                constrForceField{i}.roi.vec=fvec(checkVector,:);
            end

            dPix=50;
            min_x=min(constrForceField{i}.roi.pos(:,1))-dPix;
            max_x=max(constrForceField{i}.roi.pos(:,1))+dPix;
            min_y=min(constrForceField{i}.roi.pos(:,2))-dPix;
            max_y=max(constrForceField{i}.roi.pos(:,2))+dPix;
            maxForcePlot=1/forceField(i).par.gridSpacing*max(sqrt(forceField(i).vec(:,1).^2+forceField(i).vec(:,2).^2));
            [rows,cols]=size(currentImage);

            figure(1)
            subplot(1,2,1)
            imagesc(currentImage)
            colormap('gray')
            hold on
            quiver(forceField(i).pos(:,1),forceField(i).pos(:,2),forceField(i).vec(:,1)/maxForcePlot,forceField(i).vec(:,2)/maxForcePlot,0,'g');
            quiver(constrForceField{i}.roi.pos(:,1),constrForceField{i}.roi.pos(:,2),constrForceField{i}.roi.vec(:,1)/maxForcePlot,constrForceField{i}.roi.vec(:,2)/maxForcePlot,0,'r');
            plot(constrForceField{i}.segmRes.curve(:,1),constrForceField{i}.segmRes.curve(:,2),'k')
            plot(constrForceField{i}.segmRes.curveDilated(:,1),constrForceField{i}.segmRes.curveDilated(:,2),'r')
            % plot the holes:
            if isfield(constrForceField{i}.segmRes,'hole')
                for holeId=1:length(constrForceField{i}.segmRes.hole)
                    plot(constrForceField{i}.segmRes.hole{holeId}.curve(:,1),constrForceField{i}.segmRes.hole{holeId}.curve(:,2),'k')
                end
            end
            hold off
            title(['Constrained force field no: ',num2str(i)])
            axis equal
            xlim([max([1 min_x]) min([cols,max_x])])
            ylim([max([1 min_y]) min([rows,max_y])])        
            set(gca,'YDir','reverse')

            if strcmp(roiOK,'n') || strcmp(roiOK,'no')
                replyR=input(['Is the region of interest OK? If yes: press ENTER, else type new dilation Radius [',num2str(dilationR),']:']);
                if length(replyR)==1
                    dilationR=replyR;
                else
                    roiOK='y';
                end
            end   
        end        
        % Now that the ROI has been determined, go and find the interfaces:
        curveCellCluster=constrForceField{i}.segmRes.curveDilated;       
        
        numCellsIn=input(['How many cells do you see [',num2str(numCells),']?: ']);
        if ~isempty(numCellsIn)
            numCells=numCellsIn;
        end
        
        % Readapt the limits of the figure
        if i~=toDoList(1) && numCells>=numCellsOld && 1<=numCellsOld-1
            xlimVal=[max(min(curveInterfaceOld{1}(:,1))-50,1) min(max(curveInterfaceOld{1}(:,1))+50,cols)];
            ylimVal=[max(min(curveInterfaceOld{1}(:,2))-50,1) min(max(curveInterfaceOld{1}(:,2))+50,rows)];
        else
            xlimVal=[max([1 min_x]) min([cols,max_x])];
            ylimVal=[max([1 min_y]) min([rows,max_y])];
        end
        xlim(xlimVal)
        ylim(ylimVal)

        for j=1:numCells-1
            curveInterface=[];
            while size(curveInterface,1)<2
                display('Plot the interface')                
                if i~=toDoList(1) && numCells>=numCellsOld && j<=numCellsOld-1
                    polygonObject = impoly(gca,curveInterfaceOld{j},'Closed',false);
                    replyDummy=input('Press ENTER to continue:...');
                else
                    polygonObject = impoly(gca,'Closed',false);
                end
                curveInterface       = round(getPosition(polygonObject));
                curveInterfaceOld{j} = curveInterface;
            end
            % devide the force field into two parts:
            [cell1_mask, cell2_mask, cell1_bdr_Pix, cell2_bdr_Pix, interface_Pix, cell1_extMask, ~, ~, ~]=intersecMaskPolygon(curveCellCluster,curveInterface);
            curveCellCluster=cell2_bdr_Pix; %perimCell2_Pix;
            
            % The following lines are needed to find the inner masks:
            [rowsROI,colsROI]=size(constrForceField{i}.segmRes.mask);
            cell1_mask_large=cell1_mask;
            cell1_mask_large(rowsROI,colsROI)=0;
            cell2_mask_large=cell2_mask;
            cell2_mask_large(rowsROI,colsROI)=0;
                                   
            % store these values:
            constrForceField{i}.cell{j}.mask      = cell1_mask;
            constrForceField{i}.cell{j}.extMask   = cell1_extMask;
            constrForceField{i}.cell{j}.innerMask = constrForceField{i}.segmRes.mask &  cell1_mask_large;
            constrForceField{i}.cell{j}.center    = centerOfMass(cell1_mask);
            constrForceField{i}.cell{j}.boundary  = cell1_bdr_Pix;
            constrForceField{i}.cell{j}.interface = interface_Pix;
            constrForceField{i}.cell{j}.cellArea  = constrForceField{i}.par.pixSize_mu^2*calcCellArea(constrForceField{i}.cell{j}.innerMask,constrForceField{i}.segmRes.hole);
            % [cellArea]=um^2.

            % This will be overwritten in the next run:
            constrForceField{i}.cell{j+1}.mask      = cell2_mask;
            constrForceField{i}.cell{j+1}.extMask   = cell2_mask; % here extMask and Mask are the same
            constrForceField{i}.cell{j+1}.innerMask = constrForceField{i}.segmRes.mask &  cell2_mask_large;            
            constrForceField{i}.cell{j+1}.center    = centerOfMass(cell2_mask);
            constrForceField{i}.cell{j+1}.boundary  = cell2_bdr_Pix;
            constrForceField{i}.cell{j+1}.interface = interface_Pix;  
            constrForceField{i}.cell{j+1}.cellArea  = sum(constrForceField{i}.cell{j+1}.innerMask(:))*(constrForceField{i}.par.pixSize_mu^2);
            % [cellArea]=um^2.
            
            ROIpos=constrForceField{i}.roi.pos;
            ROIvec=constrForceField{i}.roi.vec;
            if sum(sum(abs(ROIpos-round(ROIpos))))==0
                [inposCell1,invecCell1]=findVectorFieldInMask(ROIpos,ROIvec,constrForceField{i}.cell{j}.mask);
                [inposCell2,invecCell2]=findVectorFieldInMask(ROIpos,ROIvec,constrForceField{i}.cell{j+1}.mask);
                
                % Cell1:
                constrForceField{i}.cell{j}.pos=inposCell1;
                constrForceField{i}.cell{j}.vec=invecCell1;
                
                % Cell2: This will be overwritten in the next run:
                constrForceField{i}.cell{j+1}.pos=inposCell2;
                constrForceField{i}.cell{j+1}.vec=invecCell2;
            else 
                checkVectorCell1 = inpolygon(ROIpos(:,1),ROIpos(:,2),constrForceField{i}.cell{j}.boundary(:,1)  ,constrForceField{i}.cell{j}.boundary(:,2));
                checkVectorCell2 = inpolygon(ROIpos(:,1),ROIpos(:,2),constrForceField{i}.cell{j+1}.boundary(:,1),constrForceField{i}.cell{j+1}.boundary(:,2));
                
                % Cell1:
                constrForceField{i}.cell{j}.pos=ROIpos(checkVectorCell1,:);
                constrForceField{i}.cell{j}.vec=ROIvec(checkVectorCell1,:);
                
                % Cell2: This will be overwritten in the next run:
                constrForceField{i}.cell{j+1}.pos=ROIpos(checkVectorCell2,:);
                constrForceField{i}.cell{j+1}.vec=ROIvec(checkVectorCell2,:);
            end

            
            pixSize_mu = constrForceField{i}.par.pixSize_mu;
            gridSpacing= constrForceField{i}.par.gridSpacing;
            bwMask1 = constrForceField{i}.cell{j}.mask;
            bwMask2 = constrForceField{i}.cell{j+1}.mask;
            
            % here we don't interpolate to be fast when cutting out the
            % force field. If interpolation is wanted we do this at the
            % very end!
            tic;
            [sumForceVec1,method,~,~]=integrateForceField(forceField(i).pos,forceField(i).vec,bwMask1,pixSize_mu,gridSpacing);
            [sumForceVec2,~,~,~]     =integrateForceField(forceField(i).pos,forceField(i).vec,bwMask2,pixSize_mu,gridSpacing);
            toc;
        
            % do statistics on the force Fields.
            % Sum up the force vectors for each cell:

%!!!!       % Here we multiply each stress with its support (gridsize^2*pixSize^2*10^(-12)) to
            % get the actual force in N. In the future, the support might
            % be triangles and then this becomes more complicated!
            % The conversion factor from Pa to nN is:           
            
            constrForceField{i}.cell{j}.stats.resForce.pos    =   constrForceField{i}.cell{j}.center;
            constrForceField{i}.cell{j}.stats.resForce.vec    = - sumForceVec1;        
            constrForceField{i}.cell{j}.stats.resForce.mag    =   sqrt(sum((constrForceField{i}.cell{j}.stats.resForce.vec).^2));
            constrForceField{i}.cell{j}.stats.method          =   method;

            constrForceField{i}.cell{j+1}.stats.resForce.pos    =   constrForceField{i}.cell{j+1}.center;
            constrForceField{i}.cell{j+1}.stats.resForce.vec    = - sumForceVec2;
            constrForceField{i}.cell{j+1}.stats.resForce.mag    =   sqrt(sum((constrForceField{i}.cell{j+1}.stats.resForce.vec).^2));
            constrForceField{i}.cell{j+1}.stats.method          =   method;

            % So far these interface can belong to more than one cell pair!
            constrForceField{i}.interface{j}.pos=interface_Pix.coord;
            
            % The change from resForce.vec = sum to resForce.vec = - sum 
            % might have introduced wrong signs in here, that has to be
            % checked again!
            % constrForceField{i}.interface{j}.stats.relErrorXYcomp=(constrForceField{i}.cell{j}.stats.resForce.vec+constrForceField{i}.cell{j+1}.stats.resForce.vec)./((abs(constrForceField{i}.cell{j}.stats.resForce.vec)+abs(constrForceField{i}.cell{j+1}.stats.resForce.vec))/2);
            % constrForceField{i}.interface{j}.stats.absError=(sqrt(sum((constrForceField{i}.cell{j}.stats.resForce.vec).^2))-sqrt(sum((constrForceField{i}.cell{j+1}.stats.resForce.vec).^2)))/((sqrt(sum((constrForceField{i}.cell{j}.stats.resForce.vec).^2))+sqrt(sum((constrForceField{i}.cell{j+1}.stats.resForce.vec).^2)))/2);
            % constrForceField{i}.interface{j}.stats.alpha=acos(dot(-constrForceField{i}.cell{j}.stats.resForce.vec,constrForceField{i}.cell{j+1}.stats.resForce.vec)/(constrForceField{i}.cell{j}.stats.resForce.mag*constrForceField{i}.cell{j+1}.stats.resForce.mag));

            % This is only needed to reconvert the interface force in to
            % the scales of the traction stress such that they are
            % comparable:
            factor_Pa2nN=constrForceField{i}.par.gridSpacing^2*constrForceField{i}.par.pixSize_mu^2*10^(-3);
            subplot(1,2,2)
            imagesc(currentImage)
            colormap('gray')
            hold on
            %plot innter boundary:
            plot(constrForceField{i}.segmRes.curve(:,1),constrForceField{i}.segmRes.curve(:,2),'k')
            %plot complete forceField:
            quiver(forceField(i).pos(:,1),forceField(i).pos(:,2),forceField(i).vec(:,1)/maxForcePlot,forceField(i).vec(:,2)/maxForcePlot,0,'g');            
            marker=['r','b','m','c','g','y','k'];
            for k=1:j+1                
                %plot cell{j}:
                quiver(constrForceField{i}.cell{k}.pos(:,1),constrForceField{i}.cell{k}.pos(:,2),constrForceField{i}.cell{k}.vec(:,1)/maxForcePlot,constrForceField{i}.cell{k}.vec(:,2)/maxForcePlot,0,marker(mod(k,7)+1));
                if j==1 % if there are only two cells plot the residual forces at the interface:
                    resForcePos=constrForceField{i}.cell{k}.interface.center;                    
                else % if there are more than two cells plot the residual forces at center of mass of each cell:
                    resForcePos=constrForceField{i}.cell{k}.stats.resForce.pos;
                end
%!!!            % convert the force back to a stress to plot it with the
                % traction stress. That is a bit dirty!
                quiver(resForcePos(1),resForcePos(2),constrForceField{i}.cell{k}.stats.resForce.vec(:,1)/(maxForcePlot*factor_Pa2nN),constrForceField{i}.cell{k}.stats.resForce.vec(:,2)/(maxForcePlot*factor_Pa2nN),0,marker(mod(k,7)+1));
                plot(resForcePos(1),resForcePos(2),['o',marker(mod(k,7)+1)])
                %plot boundaries:
                plot(constrForceField{i}.cell{k}.boundary(:,1),constrForceField{i}.cell{k}.boundary(:,2),marker(mod(k,7)+1))    
            end
            for k=1:j % there are only 'j' interfaces:
                plot(constrForceField{i}.interface{k}.pos(:,1),constrForceField{i}.interface{k}.pos(:,2),'--w')
            end
            % plot the holes:
            if isfield(constrForceField{i}.segmRes,'hole')
                for holeId=1:length(constrForceField{i}.segmRes.hole)
                    plot(constrForceField{i}.segmRes.hole{holeId}.curve(:,1),constrForceField{i}.segmRes.hole{holeId}.curve(:,2),'k')
                end
            end
            hold off
            title(['Constrained force field no: ',num2str(i)])
            axis equal
            xlim([max([1 min_x]) min([cols,max_x])])
            ylim([max([1 min_y]) min([rows,max_y])])        
            set(gca,'YDir','reverse')
            % save intermediate results:
            % save(path_cellCellForces, 'constrForceField','-v7.3');
            
            if j<numCells-1
                if i~=toDoList(1) && numCells>=numCellsOld && (j+1)<=numCellsOld-1
                    xlimVal=[max(min(curveInterfaceOld{j+1}(:,1))-50,1) min(max(curveInterfaceOld{j+1}(:,1))+50,cols)];
                    ylimVal=[max(min(curveInterfaceOld{j+1}(:,2))-50,1) min(max(curveInterfaceOld{j+1}(:,2))+50,rows)];
                else
                    xlimVal=[max([1 min_x]) min([cols,max_x])];
                    ylimVal=[max([1 min_y]) min([rows,max_y])];
                end
            
                subplot(1,2,1)
                imagesc(currentImage)
                colormap('gray')
                hold on
                quiver(forceField(i).pos(:,1),forceField(i).pos(:,2),forceField(i).vec(:,1)/maxForcePlot,forceField(i).vec(:,2)/maxForcePlot,0,'g');
                quiver(constrForceField{i}.cell{j+1}.pos(:,1),constrForceField{i}.cell{j+1}.pos(:,2),constrForceField{i}.cell{j+1}.vec(:,1)/maxForcePlot,constrForceField{i}.cell{j+1}.vec(:,2)/maxForcePlot,0,'r');
                plot(constrForceField{i}.segmRes.curve(:,1),constrForceField{i}.segmRes.curve(:,2),'k')
                %plot(segmRes{i}.curveDilated(:,1),segmRes{i}.curveDilated(:,2),'g')
                plot(curveCellCluster(:,1),curveCellCluster(:,2),'r')
                % mark only the relevant interfaces:
                for k=1:j % there are only 'j' interfaces:
                    on_off_list=isOnCurve(constrForceField{i}.interface{k}.pos,curveCellCluster);
                    plotPts=constrForceField{i}.interface{k}.pos(on_off_list,:);
                    plot(plotPts(:,1),plotPts(:,2),'--w');
                end
                % plot the holes:
                if isfield(constrForceField{i}.segmRes,'hole')
                    for holeId=1:length(constrForceField{i}.segmRes.hole)
                        plot(constrForceField{i}.segmRes.hole{holeId}.curve(:,1),constrForceField{i}.segmRes.hole{holeId}.curve(:,2),'k')
                    end
                end
                hold off
                title(['Constrained force field no: ',num2str(i)])
                axis equal
                xlim(xlimVal)
                ylim(ylimVal)        
                set(gca,'YDir','reverse')
            end
        end
        
        if mod(find(i==toDoList),saveInterval)==0
            display('saving the file. This might take a while:...');
            save(path_cellCellForces, 'constrForceField','-v7.3');
        end
        % This is only needed for redrawing the interfaces
        numCellsOld=numCells;
    end
        
    reply=input('All went well? Then press Enter, else type "no": ','s');
    if strcmp(reply,'no') || strcmp(reply,'n')
        finished=false;
    else
        finished=true;
    end
end

if  doIntp==1
    display('Updating constraint force field, lean back and relax:...')
    [constrForceField]=updateConstrForceField(constrForceField,forceField,'pixIntp');
end

display('saving the final results. This might take a while:...');
save(path_cellCellForces, 'constrForceField','-v7.3');
        
return;

% This has to be transfered to a different functions that plots the
% results:

startPt=1;
endPt=length(constrForceField);

figure(3)
for i=startPt:endPt
    plot(i,constrForceField{i}.stats.relErrorXYcomp(1),'or')
    hold on
    plot(i,constrForceField{i}.stats.relErrorXYcomp(2),'og')
    title('Development of rel. error over time')
end
hold off


figure(4)
for i=startPt:endPt
    plot(i,constrForceField{i}.stats.absError,'ob')
    hold on
    title('Development of abs. error over time')
end
hold off

figure(5)
for i=startPt:endPt
    plot(i,constrForceField{i}.stats.alpha*360/(2*pi),'ob')
    hold on
    title('Development of angular deviations over time')
    ylim([0 180])
end
hold off

maximumForce=0;
figure(6)
for i=startPt:endPt
    plot(constrForceField{i}.cell{1}.stats.resForce.mag,constrForceField{i}.cell{2}.stats.resForce.mag,'ob')
    hold on
    title('residual force of cell2 over residual force of cell1')
    for j=1:length(constrForceField{i}.Cell)
        currentMax=max(constrForceField{i}.cell{j}.stats.resForce.mag);
        if currentMax>maximumForce
            maximumForce=currentMax;
        end
    end
end
plot([0 maximumForce],[0 maximumForce],'--k')
hold off

figure(7)
for i=startPt:endPt
    colorStr='ob';
    for j=1:2
        if j==2
            colorStr='or';
        end
        plot(i,constrForceField{i}.cell{j}.stats.elEnergy,colorStr)
        hold on
        title('Elastic energy invested by each cell') 
    end
end
hold off

figure(8)
for i=startPt:endPt
    colorStr='ob';
    for j=1:2
        if j==2
            colorStr='or';
        end
        plot(i,constrForceField{i}.cell{j}.stats.ratioCCForceOverElEnergy,colorStr)
        hold on
        title('Ratio of cell-cell-interaction force and elastic energy') 
    end
end
hold off


function [center]=centerOfMass(BWmask)
    numPoints=sum(sum(BWmask));
    
    [rows_nest cols_nest]=size(BWmask);
    vec_x=1:cols_nest;
    mat_x=repmat(vec_x,rows_nest,1);
    
    vec_y=(1:rows_nest)';
    mat_y=repmat(vec_y,1,cols_nest);
    
    center(1)=sum(sum(BWmask.*mat_x))/numPoints;
    center(2)=sum(sum(BWmask.*mat_y))/numPoints;
end

end


