function [ output_args ] = makeNeuriteMask( analInfo ,maskDir)
%After body estimation AND reconstruction outputs the final mask 
% input will eventually be that you will be able to choose which filo you
% would like to be included in the mask. 

[upOne,~,num] = upDirectory(maskDir,1); 
%  maskDirBody =  [upOne filesep 'bodyOnly_' num]; 
%  if ~isdir(maskDirBody) 
%      mkdir(maskDirBody) 
%  end 
 % move the old files 
%  copyfile(maskDir,maskDirBody); 
if isdir(maskDir)
 rmdir(maskDir,'s'); 
end 
 mkdir(maskDir)



for iFrame = 1:length(analInfo) 
    
    maskBody = analInfo(iFrame).masks.neuriteEdge; 
    [ny,nx] = size(maskBody); 
    % make filoMask (maybe eventually make separate function) 
    filoMask   = zeros(ny,nx); 
    filoInfo = analInfo(iFrame).filoInfo; 
    if ~isempty(filoInfo) 
        
    exitFlag  =  vertcat(filoInfo(:).Ext_exitFlag); 
    type = vertcat(filoInfo(:).bodyType); %
    lengths = vertcat(filoInfo(:).Ext_length).*0.215;  
    %filoInfoPlot =   filoInfo((exitFlag >=1 & ~isnan(exitFlag)) & (type == 0 | type ==1));  %clip filoInfo by exitFlag and type
     filoInfoPlot = filoInfo((exitFlag >=1 & ~isnan(exitFlag)) & lengths > 2*0.215); 

 

 for iFilo = 1:numel(filoInfoPlot)
     
     for i = 1:numel(filoInfoPlot)
         pixIndices = filoInfoPlot(i).('Ext_pixIndices');
         idxEnd = find(pixIndices == filoInfoPlot(i).('Ext_endpointCoordFitPix'));
         pixIndicesPlot = pixIndices(1:idxEnd);
         filoMask(pixIndicesPlot) = 1;
     end
     % 
     
 end
    end 
 neuriteMask = (filoMask | maskBody); 
 
 % for now get rid of floaters (ie I have some bugs in the reconstruction
 % reattachment- better for visualization currently (20140914)
 
 neuriteMask = logical(getLargestCC(neuriteMask)); 
 %% for now do the switch my masks for whatever is located in the segmentation package so don't have to alter 
 % hunter's original function ... 
  
 
 
 
 imwrite(neuriteMask,[maskDir filesep 'neuriteFullMask' num2str(iFrame,'%03d') '.tif'],'tif'); 
 
 
 
end 

end

