function constrForceField=perfIntMeasures(constrForceField,imageFileList)

if nargin <2 || isempty(imageFileList)
   [filename, pathname] = uigetfile({'*.TIF;*.tif;*.jpg;*.png;*.*'}, ...
       'Select the first image (of e.g. Ecad)');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   imageFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for i = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end


toDoList=[];
for frame=1:length(constrForceField)
    if isfield(constrForceField{frame},'cell') && ~isempty(constrForceField{frame}.cell)
        toDoList=horzcat(toDoList,frame);
    end
end

% hard coded parameters:
r=10;

for frame=toDoList
    I = double(imread(imageFileList{frame}));

    % use different filters to average the intensity along the edge:
    % se   = strel('disk', r);
    
    h = fspecial('disk', round(r));
    % h = fspecial('average', 2*round(r)+1);
    % h = fspecial('gaussian', 6*round(r)+1, r);
    Iflt = imfilter(I,h);    
    
    % run through all interfaces and pull out the values at the pixels on
    % the interface:    
    edge=constrForceField{frame}.network.edge;
    for edgeId=1:length(edge)
        if ~isempty(edge{edgeId})
            % read out the intensity values:
            curve=edge{edgeId}.intf_internal;
            idxCurve = sub2ind(size(Iflt),curve(:,2), curve(:,1));

            edge{edgeId}.int.val   = Iflt(idxCurve);
            edge{edgeId}.int.tot   = sum(edge{edgeId}.int.val);
            edge{edgeId}.int.avg   = sum(edge{edgeId}.int.tot)/numel(idxCurve);
            edge{edgeId}.int.h     = h;
        end
    end
    constrForceField{frame}.network.edge=edge;
end