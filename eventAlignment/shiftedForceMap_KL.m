function []=shiftedForceMap_KL(inputFileList,forceField,target_dir,frameList,displField)

% max of the color scale:
cMax=[];%3000;
cMaxZoom=cMax;
cLevel=[];
% If cMax is empty, the program will take the maximal stress value.

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

if nargin < 2 || isempty(forceField)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select force field to be used as overlay');
       %the vector field:
       fileStruct=load([pathname filename]);
       forceField=fileStruct.forceField;
end

%get the target directory:
if nargin < 3 || isempty(target_dir)
    target_dir = uigetdir('','Please select target directory');
end


nImages = numel(inputFileList);
nVecFields = length(forceField);

if nargin < 4 || isempty(frameList)
    frameList=1:min(nImages,nVecFields);
end

needToCalcShift=false;
for i=frameList
    if ~isfield(forceField(i),'posShifted') || isempty(forceField(i).posShifted)
        needToCalcShift=true;
        break;
    end
end

if needToCalcShift && (nargin<5 || isempty(displField))
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select displField to be used for shifting the forceField');
       %the vector field:
       fileStruct=load([pathname filesep filename]);
       displField=fileStruct.displField;
end

if needToCalcShift
    [forceField]=shiftForceField(forceField,displField,frameList);
end

padZeros=floor(log10(length(nImages)))+1;

for i=1:min(nImages,nVecFields)
    imgInfo=imfinfo(inputFileList{1});
    maxX(i)=getfield(imgInfo,'Width');
    maxY(i)=getfield(imgInfo,'Height');
end
theXlim=max(maxX);
theYlim=max(maxY);

% Set the same force scale for all frames!
%forceScale=0;
%for i=frameList
    % Do the scaling of the quiver plot of the stress by myself:
%    maxForceComp=max(abs(forceField(i).vec(:)));
%    currForceScale=maxForceComp/forceField(i).par.gridSpacing;
%    if currForceScale>forceScale
%        forceScale=currForceScale;
%    end
%end
%forceScale=1/3*forceScale;
%displScale=2;

% calculate the length of the scale bars:
%lengthScaleBar_mu=3;
%uScaleBar_mu=lengthScaleBar_mu/displScale;
%fxScaleBar_Pa=1000;
%fyScaleBar_Pa=0;

%display('Only the first frame is plotted!')
for i=frameList
    i
    I = double(imread(inputFileList{i}));
    

    
 

  
    [rows,cols]=size(I);
    [XI,YI]=meshgrid(1:cols,1:rows);
    % for the magnitude:
    %Mblue = griddata(forceField(i).posShifted(:,1),forceField(i).posShifted(:,2),sqrt(sum(forceField(i).vec.^2,2)),XI,YI,'cubic');
    % for x-component:
    Mbluex = griddata(forceField(i).posShifted(:,1),forceField(i).posShifted(:,2),forceField(i).vec(:,1),XI,YI,'cubic');
    % for y-component:
    Mbluey = griddata(forceField(i).posShifted(:,1),forceField(i).posShifted(:,2),forceField(i).vec(:,2),XI,YI,'cubic');
    %Mblue=sqrt(Mbluex.^2+Mbluey.^2);
    % remove NaNs:
    %Mblue(isnan(Mblue))=0;
    Mbluex(isnan(Mbluex))=0;
    Mbluey(isnan(Mbluey))=0;
    
    offset=1000;
    Mbluex = Mbluex + offset;
    Mbluey = Mbluey + offset;
    
    save Bblue;
    
    %imwrite(uint16(Mblue),[target_dir,filesep,'Force_Magnitude',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiff','Compression','none') 
    imwrite(uint16(Mbluex),[target_dir,filesep,'Force_X',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiff','Compression','none') 
    imwrite(uint16(Mbluey),[target_dir,filesep,'Force_Y',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiff','Compression','none') 
end