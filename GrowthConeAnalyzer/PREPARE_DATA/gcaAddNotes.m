function [ output_args ] = gcaAddNotes(MD,plotCropped)
%gcaAddNotes: This is a simple utility function mainly for me as I need to
% parse through the data manually. It allows one to add a flag to exclude
% (which helps when making the groupLists) and make
% comments about the structure and morphology/segmentation. (Might want to
% eventually directly incorporate into movie data itself.
%%

if plotCropped == true
    plotCroppedRegion(upDirectory(MD.outputDirectory_,1));
end

s = {'LastStepChecked','Include_Currently'...
    ,'Specific_Questions','Task','Morphology','Crop','VeilStem','Filopodia','ExcludeAtStep'};


reply = listdlg('ListString',s);


% look for notes
if exist([MD.outputDirectory_ filesep 'notes.mat'],'file')~=0;
    load([MD.outputDirectory_ filesep 'notes.mat']);
else
    % initiate the structure
    notes = struct('Morphology',{''},'Crop',{''},'VeilStem',{''},'Filopodia',{''},'LastStepChecked',{''},'Include_Currently',{''},...
        'Specific_Questions',{''},'Task',{''},'ExcludeAtStep',{''});
end



toDocument = s(reply);
% put last step checked first 
% put Last Step Checked First 


for i = 1:numel(toDocument)
    % if field isn't empty show the previous comments
    
    if strcmpi(toDocument{i},'LastStepChecked');
        
        checkStepFlag  = questdlg('CheckAStep?');
        
        
       stepsFolder{1} = [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct' ...
                        filesep 'I_neurite_orientation' filesep 'Channel_1' ...
                        filesep 'CandSeeds']; 
       stepsFolder{2} = [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct'...
                        filesep 'II_neurite_orientation_refinements' filesep 'Channel_1' ];
       
       stepsFolder{3} =  [MD.outputDirectory_ filesep 'masks' filesep 'masks_for_channel_1' filesep 'Overlays'];
        
       %stepsFolder{4} =  [MD.outputDirectory_ filesep 'VisualizationOverlays/WholeNeurite/VeilWindows_SignalProcCheck/Channel_1']; 
       stepsFolder{4} = [MD.outputDirectory_ filesep 'protrusion_samples_ProtrusionBased_windSize_5ReInit61' filesep 'EdgeVelocityQuantification' ... 
           filesep 'VeilTrackCheck_ProtrusionBased_WindSize5' filesep 'Channel_1' ]; 
       
        
       % test if exist 
       idxExist = cellfun(@(x) exist(x,'file'),stepsFolder);
       idxExist = idxExist>0;
       if idxExist(2) == 1 
           load([stepsFolder{2} filesep 'backboneInfoFix.mat']); 
           if length(backboneInfo)<MD.nFrames_-1
           idxExist(2) = 0 ;
           end 
           
       end 
       stepNames{1} = 'StepI'; 
       stepNames{2}= 'StepII'; 
       stepNames{3} = 'StepIII'; 
       stepNames{4} = 'StepIV'; 
       
        while strcmpi(checkStepFlag,'yes');
            
            
            
            toCheck = stepNames(idxExist);
            replyCheck = listdlg('ListString',toCheck);
            
            
            switch replyCheck
                case 1
                    %                 cd([MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct' ...
                    %                     filesep 'I_neurite_orientation' filesep 'Channel_1' filesep 'CandSeeds']);
                    %
                    test = searchFiles('.tif',[],[MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct' ...
                        filesep 'I_neurite_orientation' filesep 'Channel_1' ...
                        filesep 'CandSeeds'],0,'all',1);
                    for iFrame = 1:size(test,1);
                        forMovie(:,:,:,iFrame) = imread(test{iFrame});
                    end
                    movieC = immovie(forMovie);
                     imSize = MD.imSize_;
                    h = setFigure(imSize(2),imSize(1),'on');
               
                    imshow(forMovie(:,:,:,1),[]);
                    implay(movieC);
                    uiwait(h);
                    
                    
                    
                case 2
                    
                   direct =   [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct'...
                        filesep 'II_neurite_orientation_refinements' filesep 'Channel_1' filesep 'AfterFix']; 
                    if exist(direct)~=0
                        
                    
                    
                    test2 = searchFiles('.tif',[],[MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct'...
                        filesep 'II_neurite_orientation_refinements' filesep 'Channel_1' filesep 'AfterFix'] ...
                        ,0,'all',1);
                    for iFrame = 1:size(test2,1);
                        forMovie2(:,:,:,iFrame) = imread(test2{iFrame});
                    end
                    movieC = immovie(forMovie2);
                      imSize = MD.imSize_;
                    h = setFigure(imSize(2),imSize(1),'on');
                  
                    imshow(forMovie2(:,:,:,1),[]);
                    
                    implay(movieC);
                    uiwait(h);
                    else 
                      msg=   msgbox('No Orientation Outliers Found');    
                      uiwait(msg); 
                    end 
                case 3 
                    direct =   [MD.outputDirectory_ filesep 'masks' filesep 'masks_for_channel_1'...
                       ]; 
                    if exist(direct)~=0
                        
                    
                    
                    test2 = searchFiles('.png',[],[MD.outputDirectory_ filesep 'masks' filesep 'masks_for_channel_1' filesep 'Overlays'] ...
                        ,0,'all',1);
                    for iFrame = 1:size(test2,1);
                        forMovie2(:,:,:,iFrame) = imread(test2{iFrame});
                    end
                    movieC = immovie(forMovie2);
                      imSize = MD.imSize_;
                    h = setFigure(imSize(2),imSize(1),'on');
                  
                    imshow(forMovie2(:,:,:,1),[]);
                    
                    implay(movieC);
                    uiwait(h);
%                     else 
%                       msg=   msgbox('No Orientation Outliers Found');    
%                       uiwait(msg); 
                    end  
                    
                case 4     
%                  direct =  [MD.outputDirectory_ filesep ...
%                      'VisualizationOverlays' filesep 'WholeNeurite' filesep 'VeilWindows_SignalProcCheck' filesep 'Channel_1']; 
                  direct= stepsFolder{4};
                  if exist(direct)~=0
                        
                    
                    
                    test2 = searchFiles('.png',[],direct ...
                        ,0,'all',1);
                    for iFrame = 1:size(test2,1);
                        forMovie2(:,:,:,iFrame) = imread(test2{iFrame});
                    end
                    movieC = immovie(forMovie2);
                      imSize = MD.imSize_;
                    h = setFigure(imSize(2),imSize(1),'on');
                  
                    imshow(forMovie2(:,:,:,1),[]);
                    
                    implay(movieC);
                    uiwait(h);
%                     else 
%                       msg=   msgbox('No Orientation Outliers Found');    
%                       uiwait(msg); 
            end 
                    
                    
            end
            
            checkStepFlag  = questdlg('CheckAStep?');
            
        end
    end
    
    
    
    
    default = notes.(toDocument{i});
            
    
    comment = inputdlg(['Comment for ' MD.outputDirectory_ ' ' toDocument{i}, '  Last Step Checked ' notes.LastStepChecked],'Comments',1,{default});
    
    notes.(toDocument{i}) = comment{1};
    save([MD.outputDirectory_ filesep 'notes.mat'],'notes');
    
    
end

