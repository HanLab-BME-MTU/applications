function [res, theta, nms, filterBank] = longsteerableDetector(img, M, sigma)

filterBank=[];
nms=[];
res=[];
theta=[];

if M==4
    resThisSigmaStack = zeros(size(img,1),size(img,2),90);
    thetaThisSigma = zeros(size(img,1),size(img,2));
  
    basicMaskSize = round(sigma*24);
     
    twoDGaussian = fspecial('gaussian',basicMaskSize*2+1,sigma);
    oneDGaussian = twoDGaussian(:,basicMaskSize);
    
    oneDLOG = -imfilter(imfilter(oneDGaussian,[1;-1],'same','replicate'),[1;-1],'same','replicate');
    
    twoDGaussian = fspecial('gaussian', basicMaskSize*2+1,5*sigma);
    oneDGaussian = twoDGaussian(basicMaskSize,:);
    Gaussian_smooth= repmat(oneDGaussian, basicMaskSize*2+1,1);
    
    twoDLineLOG = zeros(basicMaskSize*4+1,basicMaskSize*4+1);
    twoDLineLOG(basicMaskSize*1+1:basicMaskSize*3+1,basicMaskSize*1+1:basicMaskSize*3+1) ...
        = repmat(oneDLOG,1, basicMaskSize*2+1).*Gaussian_smooth;
    
    
    %          twoDLineLOG = twoDLineLOG./(sum(sum(twoDLineLOG(twoDLineLOG>0))));
    figure;imagesc(twoDLineLOG);title(['sum: ',num2str( sum(sum(twoDLineLOG(twoDLineLOG>0))))]);
    
    
    for iA = 1:90
        filter = imrotate(twoDLineLOG,iA*2,'bilinear','crop');
        %               figure(1);imagesc(filter);title(['sum: ',num2str(sum(sum(filter)))]);
        
        resThisSigmaStack(:,:,iA) = imfilter(img,filter,'same','replicate');
    end
    
    resThisSigmaStack = resThisSigmaStack.*(sigma.^(1.5));
    
    [resThisSigma,thetaThisSigma] = max(resThisSigmaStack,[],3);
    thetaThisSigma = thetaThisSigma*2;
    
    %     resAllSigma(:,:,iS) = resThisSigma;
    %     thetaAllSigma(:,:,iS) = thetaThisSigma;
    %
    
    %     figure;
    %     imagescc([resThisSigma thetaThisSigma*mean2(resThisSigma)/90]);
    %     %
    
    res = resThisSigma;
    theta = pi/2 - thetaThisSigma*pi/180;
    
end
