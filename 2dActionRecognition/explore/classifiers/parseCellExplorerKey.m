% parse key of the form 14-May-2017_m634_s13_t120_x1728_y1089_t130
function [annotationInfo] = parseCellExplorerKey(key,exp,annotation)
parsedStrings = strsplit(key,'_');

annotationInfo.exp = exp;
annotationInfo.annotationDate = parsedStrings{1};
annotationInfo.cellType = parsedStrings{2};
annotationInfo.location = parsedStrings{3};
stime =  parsedStrings{4};
annotationInfo.stime = str2double(stime(2:end));
etime = parsedStrings{7};
annotationInfo.etime = str2double(etime(2:end));
x = parsedStrings{5};
annotationInfo.x = str2double(x(2:end));
y = parsedStrings{6};
annotationInfo.y = str2double(y(2:end));
annotationInfo.annotation = annotation;
end



