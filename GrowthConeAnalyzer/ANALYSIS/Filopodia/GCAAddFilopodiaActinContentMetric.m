function [ filoInfo, avgVeilIntensity ] = GCAAddFilopodiaActinContentMetric(img,veilMask,filoInfo)
%
%% Calculate the veil intensity without internal filopodia as a normalization factor
% for the actin content intensity.
[ny,nx] =size(img); 
% filter the image based on the psf of the microscope.
H = fspecial('gaussian',3,0.43); % make variable
imgFilt = imfilter(img,H); % for weighted averaging
pixIndicesInt = cell(length(filoInfo),1);

% get the internal filo to subtract out from the veil normalization
% intensity
for i = 1:length(filoInfo)
    pixIndices = filoInfo(i).('Int_pixIndices');
    idxEnd = find(pixIndices == filoInfo(i).(['Int_endpointCoordFitPix']));
    x = pixIndices(1:idxEnd); 
    if ~isempty(x) % having trouble with concatenation with the empty cells - in future won't allow
        
    pixIndicesInt{i,1} = pixIndices(1:idxEnd);
    end
    %[yC,xC] = ind2sub( imgSize  ,pixIndicesPlot);
    % plot(xC,yC,'color',colorC,'Linewidth',3);
end
maskIntFilo = zeros(size(veilMask));
maskIntFilo(vertcat(pixIndicesInt{:})) = 1;
maskIntFilo = imdilate(maskIntFilo,(strel('disk',2))); 
veilMaskMinusFilo = (veilMask-maskIntFilo);
veilIntensityValues = imgFilt(logical(veilMaskMinusFilo));
avgVeilIntensity = mean(veilIntensityValues(:)); % normalization factor 
%%
sanityCheck = 1; 
if sanityCheck == 1 
    setFigure(nx,ny,'on'); 
    imshow(imgFilt,[]); 
    hold on 
    spy(~veilMaskMinusFilo,'b'); 
    text(5,5,['Mean Fluorescence Veil/Stem ' num2str(avgVeilIntensity,4) 'AU'],'color','y'); 
end 
    

% run through all the filopodia collect the intensity to the
% endpoint- calc an normInt Ext, normInt Int, normInt Total.
for iFilo = 1:length(filoInfo)
    
    % actin content of non-embedded filopodia first.
    pixIndices = filoInfo(iFilo).('Ext_pixIndices');
    idxEnd = find(pixIndices == filoInfo(iFilo).([ 'Ext_endpointCoordFitPix']));
    pixIndicesExt = pixIndices(1:idxEnd);
    % NOTE CHECK THE INTENSITIES THAT I SAVE IN THE FIRST STEP
    % Think for the data could likely have just saved these values
    % initially and it would save on comp time. (think the original
    % version is still the perpendicular).... 
    
    % Extra-veil filo 
    avgFiloIntensityExt =  mean(imgFilt(pixIndicesExt));
    normIntensityExt = (avgFiloIntensityExt/avgVeilIntensity);
    filoInfo(iFilo).Ext_IntensityNormToVeil = normIntensityExt;
    
    % Intra-veil filo 
    pixIndices = filoInfo(iFilo).('Int_pixIndices'); 
    idxEnd = find(pixIndices == filoInfo(iFilo).(['Int_endpointCoordFitPix']));
    pixIndicesInt = pixIndices(1:idxEnd); 
    
    avgFiloIntensityInt = mean(imgFilt(pixIndicesInt)); 
    normIntensityInt = (avgFiloIntensityInt/avgVeilIntensity);
    
    filoInfo(iFilo).Int_IntensityNormToVeil = normIntensityInt;  
  
    
    % Total Actin Bundle 
    avgFiloIntensityTot = mean(imgFilt([pixIndicesExt;pixIndicesInt])); 
    filoInfo(iFilo).Tot_IntensityNormToVeil = avgFiloIntensityTot/avgVeilIntensity; 
    
end



end

