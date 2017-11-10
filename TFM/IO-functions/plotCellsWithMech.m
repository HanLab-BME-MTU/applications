function []=plotCellsWithMech(inputFileList,vectorField,target_dir,fieldName)

%read in Stack of bead images:
if nargin < 1 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Image of the Stack to be overlaied');
   
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

fieldNames=[];
if nargin < 2 || isempty(vectorField)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select vector field to be used as overlay');
       %the vector field:
       fileStruct=load([pathname filesep filename]);
       fieldNames=char(fieldnames(fileStruct));
       if strmatch('forceField', fieldNames, 'exact')
           fieldName='forceField';
       elseif strmatch('displField', fieldNames, 'exact')
           fieldName='displField';
       else
           display('Neither displField nor forceField found')
           return
       end
       vectorField=getfield(fileStruct,fieldName);
elseif nargin < 4  % This is the case when vectorField is given but no fieldName.
   fieldName=input('What is the field Name, either "forceField" or "displField": ','s');
end

%get the target directory:
if nargin < 3 || isempty(target_dir)
    target_dir = uigetdir('','Please select target directory');
end

nImages = numel(inputFileList);
nVecFields = length(vectorField);

padZeros=floor(log10(length(nImages)))+1;

for i=1:min(nImages,nVecFields)
    imgInfo=imfinfo(inputFileList{1});
    maxX(i)=getfield(imgInfo,'Width');
    maxY(i)=getfield(imgInfo,'Height');
end
theXlim=max(maxX);
theYlim=max(maxY);

if isfield(vectorField,'par') && ~strcmp(fieldName,'displField')
        pixSize_mu=vectorField(1).par.pixSize_mu;
else
    pixSize_mu=input('Enter the pixelsize in microns: ');
end
    
for i=1:min(nImages,nVecFields)  
    % calculate the size of um in pixel:
    lengthScaleBar_mu=20;
    lengthScaleBar_pix=lengthScaleBar_mu/pixSize_mu;
    % extra spacing from the image edge in pixel:
    dPix=50;
           
    I = double(imread(inputFileList{i}));
    
    figure(1)
    colormap('gray');
    imagesc(I)
    hold on
    quiver(vectorField(i).pos(:,1),vectorField(i).pos(:,2),vectorField(i).vec(:,1),vectorField(i).vec(:,2))
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-20,[num2str(lengthScaleBar_mu),'\mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)      
    xlim([1 theXlim])
    ylim([1 theYlim])
    title(['Cells with ',fieldName])
    set(gca,'YDir','reverse','XTick',[],'YTick',[])
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.eps'], 'psc2');
    hold off
    
    Mred  =I/max(max(I));
    Mgreen=zeros(size(I));
    %Mblue =zeros(size(I));
    
    [rows,cols]=size(I);
    [XI,YI]=meshgrid(1:cols,1:rows);
    Mblue = griddata(vectorField(i).pos(:,1),vectorField(i).pos(:,2),sum(vectorField(i).vec.^2,2),XI,YI);
    Mblue(isnan(Mblue))=0;
    Mblue  =Mblue/max(max(Mblue));
    
    RGBmat(1:rows,1:cols,1)=abs(Mred-1);
    RGBmat(1:rows,1:cols,2)=Mgreen;
    RGBmat(1:rows,1:cols,3)=Mblue;
    
    figure(2)
    imagesc(RGBmat)
    title(['Cells with ',fieldName,'Magnitude'])
    imwrite(RGBmat,[target_dir,filesep,'Cells_with_',fieldName,'Magnitude',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'])  
end