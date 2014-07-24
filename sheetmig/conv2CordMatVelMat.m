function [tracksMatxCord,tracksMatyCord,velMatxCord,velMatyCord,velMatMag,checkVec]=conv2CordMatVelMat(tracksFinal,matyCord,minTrackLength,timeWindow,toDoList)
% INPUT:
% tracksFinal:  Is either matxCord containing the x-pos of the nuclei, or
%               it has the structure of movieInfo, the result of Khulouds
%               tracker.
% matyCord   :  if first argument has the structure of movieInfo, then the
%               second argument has to be empty. If first argument is a
%               matrix (defining the x-pos), then the second argument has
%               to be matrix too that defines the y-pos.

if isfield(tracksFinal,'tracksFeatIndxCG') && isempty(matyCord)
    % then the input is the tracksFinal from Khulouds tracker.
    [tracksInfoMat,tracksIndxMat] = convStruct2MatIgnoreMS(tracksFinal);

    % calculate the frame to frame velocities:
    % extract the number of columns in the Info Matrix, currently it is 8
    [~,InfoCol]=size(tracksInfoMat);
    [~,IndxCol]=size(tracksIndxMat);
    nInfoCol=InfoCol/IndxCol;
    if nInfoCol-round(nInfoCol)~=0
        dispay('Something went wrong')
        return
    end

    % extract the x and y coordinates, respectively:
    tracksMatxCord=tracksInfoMat(:,1:nInfoCol:end);
    tracksMatyCord=tracksInfoMat(:,2:nInfoCol:end);
elseif ~isstruct(tracksFinal) && numel(tracksFinal)==numel(matyCord)
    tracksMatxCord=tracksFinal;
    tracksMatyCord=matyCord;
else
    display('Function is not defined for the first argument is')
end

% take only tracks that have at least n consecutive tracked pos:
h = fspecial('average', [1,minTrackLength]);
MatxFltr=imfilter(tracksMatxCord, h, NaN);
MatyFltr=imfilter(tracksMatyCord, h, NaN);
MatSum=MatxFltr+MatyFltr;
checkVec=sum(~isnan(MatSum),2)>0;

% clean the matrices:
tracksMatxCord(~checkVec,:)=[];
tracksMatyCord(~checkVec,:)=[];

% calculate the x and y coordinates of the velcities:
velMatxCord=tracksMatxCord(:,(1+timeWindow):end)-tracksMatxCord(:,1:end-timeWindow);
velMatyCord=tracksMatyCord(:,(1+timeWindow):end)-tracksMatyCord(:,1:end-timeWindow);

velMatMag  = sqrt(velMatxCord.^2+velMatyCord.^2);

if nargin>4 && ~isempty(toDoList)
    %cut off bad frames:
    tracksMatxCord=tracksMatxCord(:,toDoList);
    tracksMatyCord=tracksMatyCord(:,toDoList);
    velMatxCord=velMatxCord(:,toDoList(1:end-timeWindow));
    velMatyCord=velMatyCord(:,toDoList(1:end-timeWindow));
    velMatMag=velMatMag(:,toDoList(1:end-timeWindow));
end