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
r1=10;
r2= 5;

for frame=toDoList
    I = double(imread(imageFileList{frame}));

    % use different filters to average the intensity along the edge:
    % se   = strel('disk', r);
    
    h1 = fspecial('disk', round(r1));
    h2 = fspecial('disk', round(r2));
    % h = fspecial('average', 2*round(r)+1);
    % h = fspecial('gaussian', 6*round(r)+1, r);
    Iflt1 = imfilter(I,h1);
    Iflt2 = imfilter(I,h2);    
    
    % run through all interfaces and pull out the values at the pixels on
    % the interface:    
    edge=constrForceField{frame}.network.edge;
    for edgeId=1:length(edge)
        if ~isempty(edge{edgeId})
            % read out the intensity values:
            curve=edge{edgeId}.intf_internal;
            idxCurve1 = sub2ind(size(Iflt1),curve(:,2), curve(:,1));
            idxCurve2 = sub2ind(size(Iflt2),curve(:,2), curve(:,1));

            edge{edgeId}.int.val    = Iflt1(idxCurve1);
            edge{edgeId}.int.val2   = Iflt2(idxCurve2);
            edge{edgeId}.int.tot    = sum(edge{edgeId}.int.val1);
            edge{edgeId}.int.tot2   = sum(edge{edgeId}.int.val2);
            edge{edgeId}.int.avg    = sum(edge{edgeId}.int.tot1)/numel(idxCurve1);
            edge{edgeId}.int.avg2   = sum(edge{edgeId}.int.tot2)/numel(idxCurve2);
            edge{edgeId}.int.h      = h1;
            edge{edgeId}.int.h2     = h2;
        end
    end
    constrForceField{frame}.network.edge=edge;
end