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

s = {'Morphology','Crop','VeilStem','Filopodia','LastStepChecked','Include_Currently'...
    ,'Specific_Questions','Task','ExcludeAtStep'};


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

for i = 1:numel(toDocument)
    % if field isn't empty show the previous comments
    
    if strcmpi(toDocument{i},'LastStepChecked');
        
        checkStepFlag  = questdlg('CheckAStep?');
        
        while strcmpi(checkStepFlag,'yes');
            
            
            
            toCheck = {'StepI','StepII','StepIII'};
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
                    
            end
            
            checkStepFlag  = questdlg('CheckAStep?');
            
        end
    end
    
    
    
    
    default = notes.(toDocument{i});
            
    
    comment = inputdlg(['Comment for ' MD.outputDirectory_ ' ' toDocument{i}, '  Last Step Checked ' notes.LastStepChecked],'Comments',1,{default});
    
    notes.(toDocument{i}) = comment{1};
    save([MD.outputDirectory_ filesep 'notes.mat'],'notes');
    
    
end

