function [elEnergy,constrForceField,constrDisplField,resForce,relDeviation]=cutOutForceFieldSingleCellOld(displField,forceField,circle,path_cell_images,path_result_dir,well_name,pauseSec)
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

if (nargin < 3) || isempty(circle)
    reply=input('Do you want to continue with default circles [0] or browse for the file [1]: ');
    if reply==1
        [filename_circle, pathname_circle] = uigetfile({'*.xls';'*.*'}, ...
        'Select the file (either .xls or .mat) containing information about the circles');
        circle=[pathname_circle,filesep,filename_circle];
    else
        for i=1:length(forceField)  
            % The center of the image is the center of the circle:
            circle(i).center(1)=min(forceField(i).pos(:,1))+round(1/2*(min(forceField(i).pos(:,1))+max(forceField(i).pos(:,1))));
            circle(i).center(2)=min(forceField(i).pos(:,2))+round(1/2*(min(forceField(i).pos(:,2))+max(forceField(i).pos(:,2))));
            % The radius is the 1/2 minimum of image width or height:
            circle(i).radius=floor(min(1/2*(max(forceField(i).pos(:,1))-min(forceField(i).pos(:,1))) , 1/2*(max(forceField(i).pos(:,2))-min(forceField(i).pos(:,2))) ) );
        end
    end
end

if ischar(circle)
    [~,~,~,ext]=getFilenameBody(circle);
    if strcmp(ext,'.xls')
        circle_xls=xlsread(circle);
        circle = convertCircle_xls_to_mat(circle_xls);
    elseif strcmp(ext,'.mat')
        fileStruct=load(circle);
        circle=fileStruct.circle;
    else
        display('Wrong file extension, expect either .xls or .mat.')
        display('Nothing has been dones.')
        return;
    end
end


%read in cell images:
if nargin < 4 || isempty(path_cell_images)
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


if (nargin < 7) || isempty(pauseSec)
    pauseSec=0;
    display(['There is a breaks of ',num2str(pauseSec),'sec between the images'])
end




n=length(circle);
padZeros=floor(log10(2*n))+1;

plotVecCircle=linspace(0,2*pi,200);
   
for i=1:n
    centerOK='n';
    radiusOK='n';    
    currentCellImage=double(imread(listCellImages{i}));
    while strcmp(centerOK,'n') || strcmp(centerOK,'no') || strcmp(radiusOK,'n') || strcmp(radiusOK,'no')
        VectorOfOnes=zeros(size(forceField(i).pos))+1;
        checkVector=(forceField(i).pos(:,1)-circle(i).center(1)).^2+(forceField(i).pos(:,2)-circle(i).center(2)).^2 < circle(i).radius^2;

        constrForceField(i).pos=forceField(i).pos(checkVector,:);
        constrForceField(i).vec=forceField(i).vec(checkVector,:);
        constrForceField(i).par=forceField(i).par;

        % Sum up the force vectors for each cell:
        resForce(i).pos=circle(i).center;
        resForce(i).vec=sum(constrForceField(i).vec,1);
        relDeviation(i)=sqrt(sum(resForce(i).vec.^2)/max(sum(constrForceField(i).vec.^2,2)));

        figure(1)
        subplot(1,2,1)
        quiver(constrForceField(i).pos(:,1),constrForceField(i).pos(:,2),constrForceField(i).vec(:,1),constrForceField(i).vec(:,2),'b');
        hold on
        plot(circle(i).center(1),circle(i).center(2),'xr')
        plot(circle(i).center(1)+circle(i).radius*cos(plotVecCircle),circle(i).center(2)+circle(i).radius*sin(plotVecCircle),'r')  
        hold off
        title(['Constrained force field no: ',num2str(i)])
        set(gca,'YDir','reverse')
        
        subplot(1,2,2)
        imagesc(currentCellImage)
        colormap('gray')
        hold on
        quiver(forceField(i).pos(:,1),forceField(i).pos(:,2),forceField(i).vec(:,1),forceField(i).vec(:,2),'b');
        plot(circle(i).center(1),circle(i).center(2),'xr')
        plot(circle(i).center(1)+circle(i).radius*cos(plotVecCircle),circle(i).center(2)+circle(i).radius*sin(plotVecCircle),'r')        
        hold off
        title(['Original force field no: ',num2str(i)])
        set(gca,'YDir','reverse')
        
        if strcmp(centerOK,'n') || strcmp(centerOK,'no')
            replyCenter=input(['Is center OK? If yes: press ENTER, else type in new center [',num2str(circle(i).center),'] :']);
            if length(replyCenter)==2
                circle(i).center(1)=replyCenter(1);
                circle(i).center(2)=replyCenter(2);
            else
                centerOK='y';
            end
        end
        
        if strcmp(radiusOK,'n') || strcmp(radiusOK,'no')
            replyRadius=input(['Is radius OK? If yes: press ENTER, else type in new radius ',num2str(circle(i).radius),' :']);
            if length(replyRadius)==1
                circle(i).radius=replyRadius;
            else
                radiusOK='y';
            end
        end        
    end
    save([pathname_forceField, filesep, 'circle.mat'],'circle');
end

display('The sum of the constrained forceField should be zero.') 
display('The sum(constrforceField)/max(constrforceField)=')
display([num2str([1:n]') repmat(': ',n,1) num2str(relDeviation')])

% Interpolate the measured displacement field on the grid of the
% constrained forcefield. Then calculate the elastic energy:

for i=1:n
    iDisplCoordx = TriScatteredInterp(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1));
    iDisplCoordy = TriScatteredInterp(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,2));

    constrDisplField(i).pos      = constrForceField(i).pos;
    constrDisplField(i).par      = displField(i).par;
    constrDisplField(i).vec(:,1) = iDisplCoordx(constrForceField(i).pos(:,1),constrForceField(i).pos(:,2));
    constrDisplField(i).vec(:,2) = iDisplCoordy(constrForceField(i).pos(:,1),constrForceField(i).pos(:,2));
    
    
    figure(10)
    subplot(1,2,1)
    quiver(constrDisplField(i).pos(:,1),constrDisplField(i).pos(:,2),constrDisplField(i).vec(:,1),constrDisplField(i).vec(:,2),'b');
    title('Interpolated and constrained displacement field')
    set(gca,'YDir','reverse')
    subplot(1,2,2)
    quiver(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1),displField(i).vec(:,2),'b');
    title('Measured displacement field')
    set(gca,'YDir','reverse')
    pause(pauseSec)
    
   
    factor=forceField(i).par.gridSpacing^2*forceField(i).par.pixSize_mu^3/10^6;    
    elEnergy(i)= 1/2*sum(sum(constrDisplField(i).vec.*constrForceField(i).vec))*factor;
end
display('The elastic energy is:');
display([num2str([1:n]') repmat(': ',n,1) num2str(elEnergy')]);


%[maxForce, meanMaxForce, stdMaxForce, MaxForceOfAllFrames]=statsForceField(constrForceField,well_name);

save([pathname_forceField, filesep, 'analysisConstrSingleCell.mat'],'elEnergy', 'constrForceField','constrDisplField','resForce','relDeviation');


if nargin>4
    if ~isdir(path_result_dir)
        mkdir(path_result_dir)
    end
    save([path_result_dir,filesep, 'analysisConstrSingleCell_',well_name,'.mat'],'elEnergy', 'constrForceField','constrDisplField','resForce','relDeviation');
end

