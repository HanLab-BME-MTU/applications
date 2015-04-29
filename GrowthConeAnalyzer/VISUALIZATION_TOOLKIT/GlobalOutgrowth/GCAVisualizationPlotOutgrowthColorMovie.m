function [ output_args ] = GCAVisualizationPlotOutgrowthColorMovie(MD)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

partitionDir = [MD.outputDirectory_ filesep 'PARAMETER_EXTRACTION\Partition_Outgrowth_Trajectory'...
    ];
neuriteOutgrowthDir  = [MD.outputDirectory_ filesep 'PARAMETER_EXTRACTION\GlobalFunctional'...
    '\neurite_outgrowth_measurements']; 


load([partitionDir filesep 'globalParams.mat']); 
pauseFrames = globalParams.outgrowth.pauseFrames; 
load([neuriteOutgrowthDir filesep 'neuriteLengthOutput.mat']); 

longPathMasks = searchFiles('.tif',[],[neuriteOutgrowthDir filesep 'Output' filesep ... 
    'LongPathMasks'],0,'all',1); 

vels = horzcat(globalParams.outgrowth.groupedVelSmoothed{:}); 
%   vels(vels>10)=10;
%    vels(vels<10) = -10; 
%                         mapper=linspace(-10,10,128)';
%                         
%                         % get closest colormap index for each feature
%                         D=createDistanceMatrix(vels,mapper);
%                         [sD,idxCMap]=sort(abs(D),2);
                        
  outDir =[ MD.outputDirectory_ filesep 'PARAMETER_EXTRACTION\GlobalFunctional\'...
    'neurite_outgrowth_measurements' filesep 'CM'];
if ~isdir(outDir)
    mkdir(outDir) 
end 
    ny = MD.imSize_(1); 
    nx = MD.imSize_(2); 
    
for iFrame = 1:MD.nFrames_
    setFigure(nx,ny,'on'); 
    img =double(imread([  MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{iFrame}])); 
    % load longPath mask 
   mask = logical(imread(longPathMasks{iFrame})); 
   imshow(-img,[]);
   hold on 
   if (vels(iFrame) <0 && isempty(intersect(iFrame,pauseFrames))); 
       c = 'b'; 
   elseif (vels(iFrame) >0 && isempty(intersect(iFrame,pauseFrames))); 
       c = 'r'; 
   else 
       c = 'g'; 
   end 
       
%    cMap=jet(128);
%    c = idxCMap(iFrame,1)
   spy(mask,c); % for now I think it is most happy with a character input - I will try to fix spy... 
   
  saveas(gcf,[outDir filesep num2str(iFrame,'%03d') '.png']); 
 
   
close gcf 
    
end 


end

