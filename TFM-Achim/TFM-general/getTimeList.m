function [timePtsRel timePts timeIntervals meanDT stdDT]=getTimeList(fileList,timeField)
% Input : list of .tif images
% Output: timePtsRel: time points in seconds relative to the first frame in the list.
%         timePts   : absolut time point in seconds.

if nargin <2 || isempty(timeField)
    timeField='FileModDate';
end

timePts=zeros(length(fileList),1);

for frame=1:length(fileList)
    imageInfo=imfinfo(fileList{frame});
    
    if strcmp(timeField,'DateTime') && isfield(imageInfo,'DateTime')
        accTime=imageInfo.DateTime;

        % This is a strange format, here we have to remove ':' so Matlab can
        % read it:
        accTime(strfind(accTime,':'))=' ';

        % convert it to a numerical vector:
        accTimeNum=str2num(accTime);

        % This can be converted now by Matlab:
        timeInDays=datenum(accTimeNum);
    elseif strcmp(timeField,'FileModDate')
        timeInDays=datenum(imageInfo.FileModDate);
    else
        display(['The field : "',timeField,'" is not supported!']);
        % all zero:
        timePtsRel=zeros(length(fileList),1);
        timePts=zeros(length(fileList),1);
        timeIntervals=zeros(length(fileList),1);
        meanDT=zeros(length(fileList),1);
        stdDT=zeros(length(fileList),1);        
        return;
    end
    
    % We convert it in seconds:
    timePts(frame) =timeInDays*24*3600;    
end
timePtsRel=timePts-timePts(1);

timeIntervals=timePts(2:end)-timePts(1:end-1);

meanDT=mean(timeIntervals);
stdDT=std(timeIntervals);

display(['The mean +- std of the time interval is: ',num2str(meanDT),'+-',num2str(stdDT)]);
display(['With [min max]= [',num2str(min(timeIntervals)),' ',num2str(max(timeIntervals)),']']);