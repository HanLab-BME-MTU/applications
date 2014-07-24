function [T,T_path]=calcStageDriftUsingCorrFlow(inputFileList,target_dir)

%read in stack of flow files:
if nargin < 1 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select the first flow file to be used');
   
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

%get the target directory:
if nargin < 2 || isempty(target_dir)
    target_dir = uigetdir('','Where to save StageDriftCorrectionByCorrFlow.mat?');
end


n = numel(inputFileList);

for i=1:n
    fileStruct=load(inputFileList{i});
    flow=fileStruct.flow;
    flow(isnan(flow(:,3)),:)=[];
    flow(isnan(flow(:,4)),:)=[];
    
    %The transformation has the same form as the registration method from
    %Sylvain. Here we take simply the median of the determined flow
    %vectors. We take the median since it is less distorted by outliers.
     
    %The backward transformation to correct for the shift in the
    %y-Coordinate is given by: 
    T(i,1)=-median(flow(:,3)-flow(:,1));
    
    %The backward transformation to correct for the shift in the
    %x-Coordinate is given by: 
    T(i,2)=-median(flow(:,4)-flow(:,2));
end
T_path=[target_dir,filesep,'StageDriftCorrectionByCorrFlow.mat'];
save(T_path, 'T');

%Not working because of flowTrack01_02 and flowTrack03_04 have different
%files bodies. Thus stack contains only one file.
% for i=1:n
%     fileStruct=load(inputFileList{i});
%     
%     flowTrack=fileStruct.flowTrack;
%     flowTrack.v{1}(isnan(flowTrack.v{1}(:,1)),:)=[];
%     flowTrack.v{1}(isnan(flowTrack.v{1}(:,2)),:)=[];
%     
%     %The transformation has the same form as the registration method from
%     %Sylvain. Here we take simply the median of the determined flow
%     %vectors. We take the median since it is less distorted by outliers.
%      
%     %The backward transformation to correct for the shift in the
%     %y-Coordinate is given by: 
%     T(i,1)=-median(flowTrack.v{1}(:,2));
%     
%     %The backward transformation to correct for the shift in the
%     %x-Coordinate is given by: 
%     T(i,2)=-median(flowTrack.v{1}(:,1));
% end
