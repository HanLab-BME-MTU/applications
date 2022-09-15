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
nz=meta.getPixelsSizeZ(0).getValue();
nc=meta.getChannelCount(0);
if nz>1 %(ip.Results.zStack==true)
    channels=reshape(raw,[],nc)';
    [r,c] = size(channels);
    conPlane=cell(1,c);
    
    for i=1:r %per channel
        inCol=channels(i,:);
        stackCol=reshape(inCol,nz,[]);
        [w, h] = size(stackCol{1,1});
        [n,m]=size(stackCol); % n is the number of layers per stack.
        meanCol=cell(1,m);
        for j=1:m % per frame
            temp = zeros(w,h,n);
            for k=1:n
                temp(:,:,n)=cell2mat(stackCol(k,j));
            end
            meanCol{j}=mean(temp,3);
        end

        conPlane(1+(i-1)*m:i*m)=meanCol; %cat(1,conPlane,
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

