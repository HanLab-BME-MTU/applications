function [constrForceField, constrDisplField]=cutOutForceFieldCells(forceField,imageFileList,target_dir)

if nargin < 1 || isempty(forceField) 
   [filename_forceField, pathname_forceField] = uigetfile({'*.mat';'*.*'}, ...
       'Select forceField.mat to be used');
       %the stage drift Transformation:
       fileStruct=load([pathname_forceField filesep filename_forceField]);
       forceField=fileStruct.forceField;
       
       target_dir=pathname_forceField;
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

finished=false;
constrForceField=[];
constrDisplField=[];
toDoList=1:length(forceField);

while finished==false    
    try
        file_path=[target_dir,filesep,'cellCellForces.mat'];
        load(file_path);
        display(['Done files are: 1...',num2str(length(constrForceField)),' of in total: 1...',num2str(length(forceField))])
        reply=input('Which fields do you want to re-analyze, Press Enter for All, or name frames as [n1 n2 ...]: ');
        if isempty(reply)
            toDoList=1:length(forceField);
            for i=toDoList
                if i<=length(constrForceField)
                    constrForceField{i}=[];
                    constrDisplField{i}=[];
                end
            end
        else
            toDoList=reply;
            for i=toDoList
                if i<=length(constrForceField)
                    constrForceField{i}=[];
                    constrDisplField{i}=[];
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
    pauseSec=0;

    %toDoList=1:51;
    for i=toDoList
        display('Plotting constraint forceField might take some time...')
        currentImage = double(imread(imageFileList{i}));
        roiOK='n';
        while strcmp(roiOK,'n') || strcmp(roiOK,'no')
            cellEdgeResults=cellPerim(imageFileList{i},dilationR,sigmaGauss,closureRadius,pauseSec);
            segmRes{i}=cellEdgeResults{1};

            checkVector = inpolygon(forceField(i).pos(:,1),forceField(i).pos(:,2),segmRes{i}.curveDilated(:,1),segmRes{i}.curveDilated(:,2));
            constrForceField{i}.roi.pos=forceField(i).pos(checkVector,:);
            constrForceField{i}.roi.vec=forceField(i).vec(checkVector,:);

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
            plot(segmRes{i}.curve(:,1),segmRes{i}.curve(:,2),'k')
            plot(segmRes{i}.curveDilated(:,1),segmRes{i}.curveDilated(:,2),'r')
            hold off
            title(['Constrained force field no: ',num2str(i)])        
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
        % Now that the ROI has been determined, go and find the interface:
        curveInterface=[];
        while size(curveInterface,1)<2
            display('Plot the interface')
            polygonObject = impoly(gca);
            curveInterface = round(getPosition(polygonObject));
        end
        % devide the force field in to two parts:            
        [~, ~, perimCell1, perimCell2,interfacePix]=intersecMaskPolygon(segmRes{i}.curveDilated,curveInterface);

        checkVector = inpolygon(forceField(i).pos(:,1),forceField(i).pos(:,2),perimCell1(:,1),perimCell1(:,2));
        constrForceField{i}.Cell{1}.pos=forceField(i).pos(checkVector,:);
        constrForceField{i}.Cell{1}.vec=forceField(i).vec(checkVector,:);

        checkVector = inpolygon(forceField(i).pos(:,1),forceField(i).pos(:,2),perimCell2(:,1),perimCell2(:,2));
        constrForceField{i}.Cell{2}.pos=forceField(i).pos(checkVector,:);
        constrForceField{i}.Cell{2}.vec=forceField(i).vec(checkVector,:);

        % store other values:   
        constrForceField{i}.Cell{1}.perim=perimCell1;
        constrForceField{i}.Cell{1}.interface=interface;
        constrForceField{i}.Cell{1}.interfacePix=interfacePix;  

        constrForceField{i}.Cell{2}.perim=perimCell2;
        constrForceField{i}.Cell{2}.interface=interface;
        constrForceField{i}.Cell{2}.interfacePix=interfacePix;

        constrForceField{i}.segmRes=segmRes{i};    
        constrForceField{i}.par=forceField(i).par;

        % do statistics on the force Fields.
        % Sum up the force vectors for each cell:
        constrForceField{i}.Cell{1}.stats.resForce.pos = constrForceField{i}.Cell{1}.interfacePix.center;
%!!!!   % To do: here we have to multiply each stress with its support (gridsize^2*pixSize^2*10^(-12)) to
        % get the actual force, so far the resForce is a stress.
        constrForceField{i}.Cell{1}.stats.resForce.vec = sum(constrForceField{i}.Cell{1}.vec,1);        
        constrForceField{i}.Cell{1}.stats.absForce     = sqrt(sum((constrForceField{i}.Cell{1}.stats.resForce.vec).^2));

        constrForceField{i}.Cell{2}.stats.resForce.pos = constrForceField{i}.Cell{2}.interfacePix.center;
        constrForceField{i}.Cell{2}.stats.resForce.vec = sum(constrForceField{i}.Cell{2}.vec,1);
        constrForceField{i}.Cell{2}.stats.absForce     = sqrt(sum((constrForceField{i}.Cell{2}.stats.resForce.vec).^2));

        % Actually these are stats specific for each interface, in the future
        % we might have more than one interfact.
        constrForceField{i}.stats.relErrorXYcomp=(constrForceField{i}.Cell{1}.stats.resForce.vec+constrForceField{i}.Cell{2}.stats.resForce.vec)./((abs(constrForceField{i}.Cell{1}.stats.resForce.vec)+abs(constrForceField{i}.Cell{2}.stats.resForce.vec))/2);
        constrForceField{i}.stats.absError=(sqrt(sum((constrForceField{i}.Cell{1}.stats.resForce.vec).^2))-sqrt(sum((constrForceField{i}.Cell{2}.stats.resForce.vec).^2)))/((sqrt(sum((constrForceField{i}.Cell{1}.stats.resForce.vec).^2))+sqrt(sum((constrForceField{i}.Cell{2}.stats.resForce.vec).^2)))/2);
        constrForceField{i}.stats.alpha=acos(dot(-constrForceField{i}.Cell{1}.stats.resForce.vec,constrForceField{i}.Cell{2}.stats.resForce.vec)/(constrForceField{i}.Cell{1}.stats.absForce*constrForceField{i}.Cell{2}.stats.absForce));

        subplot(1,2,2)
        imagesc(currentImage)
        colormap('gray')
        hold on
        %plot complete forceField:
        quiver(forceField(i).pos(:,1),forceField(i).pos(:,2),forceField(i).vec(:,1)/maxForcePlot,forceField(i).vec(:,2)/maxForcePlot,0,'g');
        %plot Cell{1}:
        quiver(constrForceField{i}.Cell{1}.pos(:,1),constrForceField{i}.Cell{1}.pos(:,2),constrForceField{i}.Cell{1}.vec(:,1)/maxForcePlot,constrForceField{i}.Cell{1}.vec(:,2)/maxForcePlot,0,'r');
        quiver(constrForceField{i}.Cell{1}.stats.resForce.pos(:,1),constrForceField{i}.Cell{1}.stats.resForce.pos(:,2),constrForceField{i}.Cell{1}.stats.resForce.vec(:,1)/maxForcePlot,constrForceField{i}.Cell{1}.stats.resForce.vec(:,2)/maxForcePlot,0,'r');
        %plot Cell{2}:
        quiver(constrForceField{i}.Cell{2}.pos(:,1),constrForceField{i}.Cell{2}.pos(:,2),constrForceField{i}.Cell{2}.vec(:,1)/maxForcePlot,constrForceField{i}.Cell{2}.vec(:,2)/maxForcePlot,0,'b');
        quiver(constrForceField{i}.Cell{2}.stats.resForce.pos(:,1),constrForceField{i}.Cell{2}.stats.resForce.pos(:,2),constrForceField{i}.Cell{2}.stats.resForce.vec(:,1)/maxForcePlot,constrForceField{i}.Cell{2}.stats.resForce.vec(:,2)/maxForcePlot,0,'b');
        %plot boundaries:
        plot(constrForceField{i}.Cell{1}.perim(:,1),constrForceField{i}.Cell{1}.perim(:,2),'r')    
        plot(constrForceField{i}.Cell{2}.perim(:,1),constrForceField{i}.Cell{2}.perim(:,2),'b')
        plot(constrForceField{i}.segmRes.curve(:,1),constrForceField{i}.segmRes.curve(:,2),'k')
        plot(constrForceField{i}.Cell{1}.interface(:,1),constrForceField{i}.Cell{1}.interface(:,2),'m')
        hold off
        title(['Constrained force field no: ',num2str(i)])        
        xlim([max([1 min_x]) min([cols,max_x])])
        ylim([max([1 min_y]) min([rows,max_y])])        
        set(gca,'YDir','reverse')
        % save intermediate results:
        save([target_dir, filesep, 'cellCellForces.mat'], 'constrForceField', 'constrDisplField');
    end

    % Solve for the displacement field of the two forcefields and then
    % calculate the elastic energy:
    for i=toDoList
        for j=1:2
            % Start with the first cell:
            %extend the functions to a regular grid:

    %!!!    %maybe instead of griddata I should use TriScatterData, might be
            %faster:
            force_x=@(x,y) myGriddata(constrForceField{i}.Cell{j}.pos(:,1),constrForceField{i}.Cell{j}.pos(:,2),constrForceField{i}.Cell{j}.vec(:,1),x,y,'cubic');
            force_y=@(x,y) myGriddata(constrForceField{i}.Cell{j}.pos(:,1),constrForceField{i}.Cell{j}.pos(:,2),constrForceField{i}.Cell{j}.vec(:,2),x,y,'cubic');

            % Here we need to calculte the displacements only at the positions where
            % the force field doesn't vanish. For the elastic energy we want to 
            % calculate int[u(x)*f(x)] and the integrand is zero if f=0.

            xmin=min(min(constrForceField{i}.Cell{j}.pos(:,1)));
            xmax=max(max(constrForceField{i}.Cell{j}.pos(:,1)));
            ymin=min(min(constrForceField{i}.Cell{j}.pos(:,2)));
            ymax=max(max(constrForceField{i}.Cell{j}.pos(:,2)));

            [ux uy]=fwdSolution(constrForceField{i}.Cell{j}.pos(:,1),constrForceField{i}.Cell{j}.pos(:,2),constrForceField{i}.par.yModu_Pa,xmin,xmax,ymin,ymax,force_x,force_y,'fft');
            constrDisplField{i}.Cell{j}.pos=constrForceField{i}.Cell{j}.pos;
            constrDisplField{i}.Cell{j}.vec=[ux uy];

            factor=forceField(i).par.gridSpacing^2*forceField(i).par.pixSize_mu^3/10^6;        
            constrForceField{i}.Cell{j}.stats.elEnergy= 1/2*sum(sum(constrDisplField{i}.Cell{j}.vec.*constrForceField{i}.Cell{j}.vec))*factor;
            constrForceField{i}.Cell{j}.stats.ratioCCForceOverElEnergy=constrForceField{i}.Cell{j}.stats.absForce/constrForceField{i}.Cell{j}.stats.elEnergy;
            clear('ux','uy')
        end
        % save intermediate results:
        save([target_dir, filesep, 'cellCellForces.mat'], 'constrForceField', 'constrDisplField');
    end
        
    reply=input('All went well? Then press Enter, else type "no": ','s');
    if strcmp(reply,'no') || strcmp(reply,'n')
        finished=false;
    else
        finished=true;
    end
end
        


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
    plot(constrForceField{i}.Cell{1}.stats.absForce,constrForceField{i}.Cell{2}.stats.absForce,'ob')
    hold on
    title('residual force of cell2 over residual force of cell1')
    for j=1:length(constrForceField{i}.Cell)
        currentMax=max(constrForceField{i}.Cell{j}.stats.absForce);
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
        plot(i,constrForceField{i}.Cell{j}.stats.elEnergy,colorStr)
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
        plot(i,constrForceField{i}.Cell{j}.stats.ratioCCForceOverElEnergy,colorStr)
        hold on
        title('Ratio of cell-cell-interaction force and elastic energy') 
    end
end
hold off


function [Zout]=myGriddata(Xin,Yin,Zin,Xout,Yout,method)
    Zout=griddata(Xin, Yin, Zin, Xout, Yout, method);
    nanMat=isnan(Zout);
    Zout(nanMat)=0;
end

end


