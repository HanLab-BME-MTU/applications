function []=compare_BEM_FTTC_forceFields(inputFileList,ffBEM,ffFTTC,frameList,displField)

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

if nargin < 2 || isempty(ffBEM)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select force field to be used as overlay');
       %the vector field:
       fileStruct=load([pathname filesep filename]);
       ffBEM=fileStruct.forceField;
end

if nargin < 3 || isempty(ffFTTC)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select force field to be used as overlay');
       %the vector field:
       fileStruct=load([pathname filesep filename]);
       ffFTTC=fileStruct.forceField;
end

nImages = numel(inputFileList);
nVecFields = length(ffBEM);

if nargin < 4 || isempty(frameList)
    frameList=1:min(nImages,nVecFields);
end

needToCalcShift=false;
for i=frameList
    if isempty(ffBEM(i).posShifted) || isempty(ffFTTC(i).posShifted)
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
    [ffBEM] =shiftForceField( ffBEM,displField,frameList);
    [ffFTTC]=shiftForceField(ffFTTC,displField,frameList);
end

for i=1:min(nImages,nVecFields)
    imgInfo=imfinfo(inputFileList{1});
    maxX(i)=getfield(imgInfo,'Width');
    maxY(i)=getfield(imgInfo,'Height');
end
theXlim=max(maxX);
theYlim=max(maxY);


display('Only the first frame is plotted!')
for i=frameList
    % in order to avoid that the force field will be shifted twice, we
    % check if it has been shifted before:
    
    pixSize_mu=ffBEM(i).par.pixSize_mu;
    
    % calculate the size of um in pixel:
    lengthScaleBar_mu=20;
    lengthScaleBar_pix=lengthScaleBar_mu/pixSize_mu;
    % extra spacing from the image edge in pixel:
    dPix=50;
           
    I = double(imread(inputFileList{i}));  
    
    
    
    % calculate the magnitude of the force Fields:
     ffBEM(i).mag=sqrt(sum( ffBEM(i).vec.^2,2));
    ffFTTC(i).mag=sqrt(sum(ffFTTC(i).vec.^2,2));
    
    maxForceMag=max(max(ffBEM(i).mag),max(ffFTTC(i).mag));
    forceScale=0.1*maxForceMag;
    
   
    % This is the direct comparison of the two force fields:
    figure(1)
    colormap('gray');
    imagesc(I)
    hold on
    quiver(ffBEM(i).posShifted(:,1),ffBEM(i).posShifted(:,2),ffBEM(i).vec(:,1)/forceScale,ffBEM(i).vec(:,2)/forceScale,0,'b')
    quiver(ffFTTC(i).posShifted(:,1),ffFTTC(i).posShifted(:,2),ffFTTC(i).vec(:,1)/forceScale,ffFTTC(i).vec(:,2)/forceScale,0,'r')
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-20,[num2str(lengthScaleBar_mu),'\mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)      
    xlim([1 theXlim])
    ylim([1 theYlim])
    title('Comparison of the two force fields')
    set(gca,'YDir','reverse')%,'XTick',[],'YTick',[])
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.eps'], 'psc2');
    hold off
    
    cutoff=0;
     ffBEMHist= ffBEM(i).mag( ffBEM(i).mag>cutoff);
    ffFTTCHist=ffFTTC(i).mag(ffFTTC(i).mag>cutoff);
    figure(2)
    subplot(1,2,1)
    hist(ffBEMHist,1000)
    title('BEM: distribution of force magnitude')
    xlim([1,maxForceMag]);
    subplot(1,2,2)
    hist(ffFTTCHist,1000)
    title('FTTC: distribution of force magnitude')
    xlim([1,maxForceMag]);
    display(['Mean of BEM:  ',num2str(mean( ffBEMHist))]);
    display(['Mean of FTTC: ',num2str(mean(ffFTTCHist))]);
    
    display(['Median of BEM:  ',num2str(median( ffBEMHist))]);
    display(['Median of FTTC: ',num2str(median(ffFTTCHist))]);
    
    display(['# points of BEM:  ',num2str(length( ffBEMHist))]);
    display(['# points of FTTC: ',num2str(length(ffFTTCHist))]);
    
    
    [rows,cols]=size(I);
    [XI,YI]=meshgrid(1:cols,1:rows);
    imag = griddata(ffBEM(i).posShifted(:,1),ffBEM(i).posShifted(:,2),ffBEM(i).mag,XI,YI,'cubic');
    imag(isnan(imag))=0;
    ffBEM(i).imag=imag;
    
    [rows,cols]=size(I);
    [XI,YI]=meshgrid(1:cols,1:rows);
    imag = griddata(ffFTTC(i).posShifted(:,1),ffFTTC(i).posShifted(:,2),ffFTTC(i).mag,XI,YI,'cubic');
    imag(isnan(imag))=0;
    ffFTTC(i).imag=imag;
    
    figure(3)
    subplot(1,2,1)
    colormap('jet')
    imagesc(ffBEM(i).imag)
    title('BEM: force magnitude')
    caxis([0 maxForceMag])
    subplot(1,2,2)
    colormap('jet')
    imagesc(ffFTTC(i).imag)
    caxis([0 maxForceMag])
    title('FTTC: force magnitude')
end    
