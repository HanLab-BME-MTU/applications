function [fieldGeom]=createFieldGeom(leftUpperCorner,lowerRightCorner,edgeErrode,target_path)

if nargin>2 && ~isempty(edgeErrode)
    leftUpperCorner=leftUpperCorner+edgeErrode;
    lowerRightCorner=lowerRightCorner-edgeErrode;
end

fieldGeom{1}.MP=leftUpperCorner;
fieldGeom{1}.label='F1';

fieldGeom{1}.bndX(1,1)=leftUpperCorner(1);
fieldGeom{1}.bndX(2,1)=leftUpperCorner(1);
fieldGeom{1}.bndX(3,1)=lowerRightCorner(1);
fieldGeom{1}.bndX(4,1)=lowerRightCorner(1);
fieldGeom{1}.bndX(5,1)=leftUpperCorner(1);

fieldGeom{1}.bndY(1,1)=leftUpperCorner(2);
fieldGeom{1}.bndY(2,1)=lowerRightCorner(2);
fieldGeom{1}.bndY(3,1)=lowerRightCorner(2);
fieldGeom{1}.bndY(4,1)=leftUpperCorner(2);
fieldGeom{1}.bndY(5,1)=leftUpperCorner(2);

fieldGeom{1}.gridDx=[];
fieldGeom{1}.gridDy=[];
fieldGeom{1}.gridPts=[];

if (nargin<4) || isempty(target_path)
    save('fieldGeom.mat','fieldGeom');
    display('The fieldGeom.mat file has been saved to the current directory!')
elseif isdir(target_path)
    save([target_path,filesep,'fieldGeom.mat'],'fieldGeom');
else
    display('Nothing has been done, something went wrong')
end