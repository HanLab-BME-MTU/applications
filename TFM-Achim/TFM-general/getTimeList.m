function [timePtsRel timePts timeIntervals meanDT stdDT]=getTimeList(fileList,timeField,batch)
% Input : list of .tif images
% Output: timePtsRel: time points in seconds relative to the first frame in the list.
%         timePts   : absolut time point in seconds.
doBoth=0;
if nargin <2 || isempty(timeField)
    doBoth=1;
end
if nargin<3 || isempty(batch)
    batch=0;
end

timePts=zeros(length(fileList),1);


imageInfo=imfinfo(fileList{1});
if (doBoth && (~isfield(imageInfo,'DateTime') || ~isfield(imageInfo,'FileModDate'))) ||...
    ~isfield(imageInfo,timeField)
    timePtsRel   =NaN*timePts;
    timePts      =NaN*timePts;
    timeIntervals=NaN*timePts;
    meanDT       =NaN;
    stdDT        =NaN;
    display('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    display('!!! could not find all necessary timefields!!!')
    display('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    return;
end
    
foundDateTime=0;
foundFileModDate=0;


for frame=1:length(fileList)
    imageInfo=imfinfo(fileList{frame});
    
    if (doBoth && isfield(imageInfo,'DateTime')) ||...
            (strcmp(timeField,'DateTime') && isfield(imageInfo,'DateTime'))
        accTime=imageInfo.DateTime;
        
        % This can be converted now by Matlab:
        try
            timeInDays_DateTime=datenum(accTime,'yyyymmdd HH:MM:SS.FFF');
        catch
            timeInDays_DateTime=datenum(accTime,'yyyy:mm:dd HH:MM:SS');
        end

        % The date could be recovered by e.g.:
        % datestr(timeInDays_DateTime,'yyyy mm dd HH:MM:SS.FFF')
        foundDateTime=1;
    end
    if (doBoth && isfield(imageInfo,'FileModDate')) ||...
            (strcmp(timeField,'FileModDate') && isfield(imageInfo,'FileModDate'))
        timeInDays_FileModDate=datenum(imageInfo.FileModDate);
        foundFileModDate=1;
    end
    
    if ~foundFileModDate && ~foundDateTime
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
    if foundFileModDate
        timePts_FileModDate(frame) =timeInDays_FileModDate*24*3600;
    end
    if foundDateTime
        timePts_DateTime(frame) =timeInDays_DateTime*24*3600;
    end
end
if foundFileModDate
    timePtsRel_FileModDate=timePts_FileModDate-timePts_FileModDate(1);
    timeIntervals_FileModDate=timePts_FileModDate(2:end)-timePts_FileModDate(1:end-1);
    meanDT_FileModDate=mean(timeIntervals_FileModDate);
    stdDT_FileModDate=std(timeIntervals_FileModDate);
end
if foundDateTime
    timePtsRel_DateTime=timePts_DateTime-timePts_DateTime(1);
    timeIntervals_DateTime=timePts_DateTime(2:end)-timePts_DateTime(1:end-1);
    meanDT_DateTime=mean(timeIntervals_DateTime);
    stdDT_DateTime=std(timeIntervals_DateTime);
end

if foundFileModDate && foundDateTime
    % check that the values are the same! The check is positive if the
    % root mean squared deviation between the two time series is less then
    % the sum of both standard deviations.
    if sqrt(mean((timeIntervals_FileModDate-timeIntervals_DateTime).^2))<(stdDT_FileModDate+stdDT_DateTime)
        timePtsRel    = timePtsRel_DateTime;
        timePts       = timePts_DateTime;
        timeIntervals = timeIntervals_DateTime;
        meanDT        = meanDT_DateTime;
        stdDT         = stdDT_DateTime;
    else
        display('**************************************')
        display(['meanDT_DateTime:= ',num2str(meanDT_DateTime),'+-',num2str(stdDT_DateTime),'   meanDT_FileModDate:= ',num2str(meanDT_FileModDate),'+-',num2str(stdDT_FileModDate)]);
        display('**************************************')
        display('* The extracted values do not agree! *') 
        display('**************************************')
        if ~batch
            input('Should we go on with the DateTime values? Hit enter, else abort: ');
        end
        
        timePtsRel    = timePtsRel_DateTime;
        timePts       = timePts_DateTime;
        timeIntervals = timeIntervals_DateTime;
        meanDT        = meanDT_DateTime;
        stdDT         = stdDT_DateTime;        
    end
elseif foundFileModDate && ~foundDateTime
        timePtsRel    = timePtsRel_FileModDate;
        timePts       = timePts_FileModDate;
        timeIntervals = timeIntervals_FileModDate;
        meanDT        = meanDT_FileModDate;
        stdDT         = stdDT_FileModDate;
else
        timePtsRel    = timePtsRel_DateTime;
        timePts       = timePts_DateTime;
        timeIntervals = timeIntervals_DateTime;
        meanDT        = meanDT_DateTime;
        stdDT         = stdDT_DateTime;
end
    


display(['The mean +- std of the time interval is: ',num2str(meanDT),'+-',num2str(stdDT)]);
display(['With [min max]= [',num2str(min(timeIntervals)),' ',num2str(max(timeIntervals)),']']);