function ND2TIFF(varargin)

ip = inputParser;
ip.addRequired('inputPath', @ischar);
ip.addRequired('outputPath', @ischar);
ip.addParameter('zStack', false , @islogical);
ip.parse(varargin{:});
%'H:\han-research\data testing\02\nk _live_02.nd2
ret = bfopen(ip.Results.inputPath);
raw=ret{1,1}(:,1);
meta=ret{1,4};
conPlane=[];
nz=meta.getPixelsSizeZ(0).getValue();
nc=meta.getChannelCount(0);
if(ip.Results.zStack==true)
    channels=reshape(raw,[],nc)';
    [r,c] = size(channels);
    
    for i=1:r
        inCol=channels(1,:);
        stackCol=reshape(inCol,nz,[]);
        [n,m]=size(stackCol);
        meanCol=[];
        for j=1:m
            temp=uint16(zeros(512,512));
            for k=1:n
                temp=temp+reshape(cell2mat(stackCol(k,j)),512,512);
                %pause(0.01);
            end
            temp={temp./nz};
            meanCol=[meanCol;temp];
        end
        conPlane=cat(1,conPlane,meanCol);
    end
    toInt = @(x) javaObject('ome.xml.model.primitives.PositiveInteger', ...
                        javaObject('java.lang.Integer', x));
    meta.setPixelsSizeZ(toInt(1),0);
    nz=1;


end
stackPlane=reshape(reshape(conPlane,[],nc)',[],1);
planes=reshape(stackPlane,1,1,nc,nz,[]);
matPlanes=cell2mat(planes);
bfsave(matPlanes,ip.Results.outputPath,'XYCZT','metadata',meta);



%function ND2TIFF(varargin)

% ip = inputParser;
% ip.addRequired('inputPath', @ischar);
% ip.addRequired('outputPath', @ischar);
% ip.addParameter('zStack', false , @islogical);
% ip.parse(varargin{:});
% %'H:\han-research\data testing\02\nk _live_02.nd2
% ret = bfopen(ip.Results.inputPath);
% meta=ret{1,4};
% nz=meta.getPixelsSizeZ(0).getValue();
% planes=reshape(ret{1,1}(:,1),1,1,4,nz,[]);
% matPlanes=cell2mat(planes);
% if(ip.Results.zStack==true)
%     meanMatPlanes=mean(matPlanes,4);
%     meanMatPlanes
%     toInt = @(x) javaObject('ome.xml.model.primitives.PositiveInteger', ...
%                         javaObject('java.lang.Integer', x));
%     meta.setPixelsSizeZ(toInt(1),0);
% end
% bfsave(matPlanes,ip.Results.outputPath,'XYCZT','metadata',meta);

