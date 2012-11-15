function [maxResp,d2X,d2Y,d2Z,maxRespScale] = multiscaleSurfaceFilter3D(imageIn,varargin)


ip = inputParser;
%ip.addRequired('imageIn',@(x)(ndims(x) == 3));
ip.addParamValue('SigmasXY',[1 2 4],@(x)(all(x>=1)));
ip.addParamValue('SigmasZ',[1 2 4],@(x)(all(x>=1)));
ip.addParamValue('WeightZ',2,@(x)(numel(x) == 1 && x > 0));
ip.addParamValue('NormResp',false,@(x)(numel(x) == 1 && islogical(x)));
ip.parse(varargin{:});
p = ip.Results;

% sigmasXY = [1.5 2 4 ];%TTTEEEEMMMPPPP!!!
% sigmasZ = sigmasXY*4/3;

nSig = numel(p.SigmasXY);

maxResp = zeros(size(imageIn));
maxRespScale = zeros(size(imageIn));
d2X = zeros(size(imageIn));
d2Y = zeros(size(imageIn));
d2Z = zeros(size(imageIn));

for j = 1:nSig
            
    [d2Xtmp,d2Ytmp,d2Ztmp] = surfaceFilterGauss3D(imageIn,[p.SigmasXY(j) p.SigmasXY(j) p.SigmasZ(j)]);            
    
    d2Xtmp(d2Xtmp<0) = 0;
    d2Ytmp(d2Ytmp<0) = 0;
    d2Ztmp(d2Ztmp<0) = 0;        
        
    
    if p.NormResp
        %Get magnitude and normalize response based on sigma to give comparable
        %responses at different scales. This isn't even the right way to normalize, but we leave it for now in case we need to reproduce old results. Someday I'll have the time to fix this and put in proper scale-normalization....hahahahahahahahah that's funny        
        sMag = sqrt((d2Xtmp * p.SigmasXY(j)) .^2 + ...
                    (d2Ytmp * p.SigmasXY(j)) .^2 + ...
                    (d2Ztmp * p.SigmasZ(j)) .^2 .* p.WeightZ);
    else
        
        sMag = sqrt(d2Xtmp .^2 + d2Ytmp .^2 + d2Ztmp .^2);
    end
    
    isBetter = sMag > maxResp;    
    maxRespScale(isBetter) = j; 
    maxResp(isBetter) = sMag(isBetter);
    d2X(isBetter) = d2Xtmp(isBetter);
    d2Y(isBetter) = d2Ytmp(isBetter);
    d2Z(isBetter) = d2Ztmp(isBetter);
                
end