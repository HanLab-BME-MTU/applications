function [ output_args ] = gcaAddNotes(MD)
%gcaAddNotes: This is a simple utility function mainly for me as I need to 
% parse through the data manually. It allows one to add a flag to exclude
% (which helps when making the groupLists) and make
% comments about the structure and morphology/segmentation. (Might want to 
% eventually directly incorporate into movie data itself. 
%%
s = {'Morphology','Crop','VeilStem','Filopodia','Include_Currently'...
   ,'Specific_Questions','Task','ExcludeAtStep'}; 


reply = listdlg('ListString',s); 

% look for notes 
if exist([MD.outputDirectory_ filesep 'notes.mat'],'file')~=0; 
    load([MD.outputDirectory_ filesep 'notes.mat']); 
else 
    % initiate the structure 
   notes = struct('Morphology',{''},'Crop',{''},'VeilStem',{''},'Filopodia',{''},'Include_Currently',{''},...
    'Specific_Questions',{''},'Task',{''},'ExcludeAtStep',{''}); 
end 

toDocument = s(reply); 
    
for i = 1:numel(toDocument)
    % if field isn't empty show the previous comments 

default = notes.(toDocument{i}); 

comment = inputdlg(['Comment for ' toDocument{i}],'Comments',1,{default});

notes.(toDocument{i}) = comment{1}; 
save([MD.outputDirectory_ filesep 'notes.mat'],'notes'); 


end

