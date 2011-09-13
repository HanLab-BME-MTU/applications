function [constrForceField,constrDisplField]=cutOutForceFieldSingleCell(displField,forceField,path_cell_images,target_dir,pauseSec) % this should be done in TFM_part_3... :path_result_dir,well_name)
display('The energy is only correct for regular square lattices!')

if (nargin < 1) || isempty(displField)
   [filename_displField, pathname_displField] = uigetfile({'*.mat';'*.*'}, ...
       'Select displField.mat to be used');
       fileStruct=load([pathname_displField filesep filename_displField]);
       displField=fileStruct.displField;
elseif ischar(displField)
       fileStruct=load(displField);
       displField=fileStruct.displField;
end

if (nargin < 2) || isempty(forceField)
   [filename_forceField, pathname_forceField] = uigetfile({'*.mat';'*.*'}, ...
       'Select forceField.mat to be used');
       fileStruct=load([pathname_forceField filesep filename_forceField]);
       forceField=fileStruct.forceField;
elseif ischar(forceField)
       [path,~,~,~]=getFilenameBody(forceField,filesep);
       pathname_forceField=path;

       fileStruct=load(forceField);
       forceField=fileStruct.forceField;
end

%read in cell images:
if nargin < 3 || isempty(path_cell_images)
   [filename, pathname] = uigetfile({'*.tif';'*.*'}, ...
       'Select first cell image to be overlaid');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end   
   listCellImages = getFileStackNames([pathname filesep filename]);
elseif ischar(path_cell_images)
   listCellImages=getFileListFromFolder(path_cell_images);
elseif iscell(path_cell_images) 
   listCellImages=path_cell_images;
end

%get the target directory:
if nargin < 4 || isempty(target_dir)
    target_dir = uigetdir('','Select target directory for singleCellForces.mat');
end

if (nargin < 5) || isempty(pauseSec)
    pauseSec=0;
    display(['There is a breaks of ',num2str(pauseSec),'sec between the images'])
end

n=length(forceField);
finished=false;
constrForceField=[];
constrDisplField=[];
toDoList=1:n;
padZeros=floor(log10(2*n))+1;

while finished==false
    try
        file_path=[target_dir,filesep,'singleCellForces.mat'];
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
            display('No "singleCellForces.mat" file found, run through analysis completely.')
        else
            display(['Unknown error: ',exception.identifier])
        end
    end
    
    for i=toDoList
        currentImage=double(imread(listCellImages{i}));
        display(['Draw the ellipse for frame ',num2str(i),' finish with double-click: '])
        
        % For plotting the vectors in the right scale:
        maxForcePlot=1/forceField(i).par.gridSpacing*max(sqrt(forceField(i).vec(:,1).^2+forceField(i).vec(:,2).^2));
        
        figure(1)
        subplot(1,2,1)
        h_im=imagesc(currentImage);        
        colormap('gray');
        title(['Original force field no: ',num2str(i)])
        hold on
        quiver(forceField(i).pos(:,1),forceField(i).pos(:,2),forceField(i).vec(:,1)/maxForcePlot,forceField(i).vec(:,2)/maxForcePlot,0,'g');
        if i==1
            % Create a default ellipse for the first frame. % The center of
            % the image is the center of the ellipse, its radius is 
            % min(100,ImageWidth/2,ImageHeight/2), so we start with a circle:
            center(1)=min(forceField(i).pos(:,1))+round(1/2*(min(forceField(i).pos(:,1))+max(forceField(i).pos(:,1))));
            center(2)=min(forceField(i).pos(:,2))+round(1/2*(min(forceField(i).pos(:,2))+max(forceField(i).pos(:,2))));
            radius=floor(min([100, 1/2*(max(forceField(i).pos(:,1))-min(forceField(i).pos(:,1))) , 1/2*(max(forceField(i).pos(:,2))-min(forceField(i).pos(:,2)))]));
            
            % The ellipse is defined by a rectangle with [xmin ymin width height]:
            constrForceField{i}.ellipse.pos=[center(1)-radius, center(2)-radius, radius, radius];
            
            % Plot the ellipse at this position:
            h_el = imellipse(gca, constrForceField{i}.ellipse.pos);
            wait(h_el);
        else
            % Else, we take the ellipse from the previous frame as default:
            h_el = imellipse(gca, constrForceField{i-1}.ellipse.pos);
            wait(h_el);
        end
        hold off
        set(gca,'YDir','reverse')
        
        % Read out the position and a mask for the drawn ellipse:
        constrForceField{i}.ellipse.pos  = getPosition(h_el);        % This generates the position rectangle.
        constrForceField{i}.ellipse.mask =  createMask(h_el,h_im);   % This generates a BW mask
        
        % This generates the curve of the ellipse using the mask:
        B=bwboundaries(constrForceField{i}.ellipse.mask,'noholes');
        constrForceField{i}.ellipse.curve(:,1)=B{1}(:,2);
        constrForceField{i}.ellipse.curve(:,2)=B{1}(:,1);        
        
        checkVector = inpolygon(forceField(i).pos(:,1),forceField(i).pos(:,2),constrForceField{i}.ellipse.curve(:,1),constrForceField{i}.ellipse.curve(:,2));
        constrForceField{i}.pos=forceField(i).pos(checkVector,:);
        constrForceField{i}.vec=forceField(i).vec(checkVector,:);
        constrForceField{i}.par=forceField(i).par;    
        % Sum up the force vectors for each cell:
        constrForceField{i}.stats.resForce.pos(1)=constrForceField{i}.ellipse.pos(1)+0.5*constrForceField{i}.ellipse.pos(3);
        constrForceField{i}.stats.resForce.pos(2)=constrForceField{i}.ellipse.pos(2)+0.5*constrForceField{i}.ellipse.pos(4);
        constrForceField{i}.stats.resForce.vec   =sum(constrForceField{i}.vec,1);
        constrForceField{i}.stats.relDeviation   =sqrt(sum(constrForceField{i}.stats.resForce.vec.^2)/max(sum(constrForceField{i}.vec.^2,2)));
        
        
        dPix=50;
        min_x=min(constrForceField{i}.pos(:,1))-dPix;
        max_x=max(constrForceField{i}.pos(:,1))+dPix;
        min_y=min(constrForceField{i}.pos(:,2))-dPix;
        max_y=max(constrForceField{i}.pos(:,2))+dPix;
        [rows,cols]=size(currentImage);

        subplot(1,2,2)
        imagesc(currentImage)
        colormap('gray')
        hold on
        quiver(forceField(i).pos(:,1),forceField(i).pos(:,2),forceField(i).vec(:,1)/maxForcePlot,forceField(i).vec(:,2)/maxForcePlot,0,'g');
        quiver(constrForceField{i}.pos(:,1),constrForceField{i}.pos(:,2),constrForceField{i}.vec(:,1)/maxForcePlot,constrForceField{i}.vec(:,2)/maxForcePlot,0,'r');
        quiver(constrForceField{i}.stats.resForce.pos(:,1),constrForceField{i}.stats.resForce.pos(:,2),constrForceField{i}.stats.resForce.vec(1)/maxForcePlot,constrForceField{i}.stats.resForce.vec(2)/maxForcePlot,0,'k');
        plot(constrForceField{i}.ellipse.curve(:,1),constrForceField{i}.ellipse.curve(:,2),'r')
        hold off
        title(['Constrained force field no: ',num2str(i)])        
        xlim([max([1 min_x]) min([cols,max_x])])
        ylim([max([1 min_y]) min([rows,max_y])])        
        set(gca,'YDir','reverse')
        
        % save intermediate results:
        save([target_dir, filesep, 'singleCellForces.mat'], 'constrForceField', 'constrDisplField');
    end
    
    display('The sum of the constrained forceField should be zero.') 
    display('The sum(constrforceField)/max(constrforceField)=')
    numFields=length(constrForceField);
    for i=1:numFields
        allDeviations(i)=constrForceField{i}.stats.relDeviation;
    end    
    display([num2str((1:numFields)') repmat(': ',numFields,1) num2str(allDeviations')])

    
    reply=input('All went well? Then press Enter, else type "no": ','s');
    if strcmp(reply,'no') || strcmp(reply,'n')
        finished=false;
    else
        finished=true;
    end
end


% Interpolate the measured displacement field on the grid of the
% constrained forcefield. Then calculate the elastic energy:

for i=1:numFields
    currentImage=double(imread(listCellImages{i}));
    
    iDisplCoordx = TriScatteredInterp(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1));
    iDisplCoordy = TriScatteredInterp(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,2));

    constrDisplField{i}.pos      = constrForceField{i}.pos;   
    constrDisplField{i}.vec(:,1) = iDisplCoordx(constrForceField{i}.pos(:,1),constrForceField{i}.pos(:,2));
    constrDisplField{i}.vec(:,2) = iDisplCoordy(constrForceField{i}.pos(:,1),constrForceField{i}.pos(:,2));
    constrDisplField{i}.par      = displField(i).par;
    
    
    figure(10)
    subplot(1,2,1)
    imagesc(currentImage);        
    colormap('gray');
    hold on
    quiver(constrDisplField{i}.pos(:,1),constrDisplField{i}.pos(:,2),constrDisplField{i}.vec(:,1),constrDisplField{i}.vec(:,2),'r');
    plot(constrForceField{i}.ellipse.curve(:,1),constrForceField{i}.ellipse.curve(:,2),'r')
    title(['Interp. and constrained displ field, frame: ',num2str(i)])
    set(gca,'YDir','reverse')
    hold off
    
    subplot(1,2,2)
    imagesc(currentImage);        
    colormap('gray');
    hold on
    quiver(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1),displField(i).vec(:,2),'b');
    plot(constrForceField{i}.ellipse.curve(:,1),constrForceField{i}.ellipse.curve(:,2),'r')
    title(['Measured displ field, frame: ',num2str(i)])
    set(gca,'YDir','reverse')
    hold off
    
    pause(pauseSec);
       
    factor=forceField(i).par.gridSpacing^2*forceField(i).par.pixSize_mu^3/10^6;    
    constrForceField{i}.stats.elEnergy= 1/2*sum(sum(constrDisplField{i}.vec.*constrForceField{i}.vec))*factor;
end
display('The elastic energy is:');
for i=1:numFields
    allelEnergy(i)=constrForceField{i}.stats.elEnergy;
end  
display([num2str((1:numFields)') repmat(': ',numFields,1) num2str(allelEnergy')]);

save([target_dir, filesep, 'singleCellForces.mat'], 'constrForceField', 'constrDisplField');

