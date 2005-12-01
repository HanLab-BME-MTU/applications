function [tRes,sigD] = adgui_calcPlotData_testSignificance(anaDat,tags,dataProperties)
%calculate whether a distance change is significant knowing the uncertainities
%of the positional estimates (0 if not, +1 if distance increases or decreases (this calculation
% is done in pix, so increase or decrease might not necessarily be reflected correctly!))

%init vars & consts
tRes=zeros(length(anaDat)-1,1);
sigD = zeros(size(tRes));
TEST_PROB=0.01;

distanceMatrix=cat(3,anaDat.distanceMatrix);
distance = squeeze(distanceMatrix(tags(1),tags(2),:));

distanceVectorMatrixN = cat(4, anaDat.distanceVectorMatrixN);
distanceVector = squeeze(distanceVectorMatrixN(tags(1),tags(2),:,:))' .* (distance*ones(1,3));



%loop through all timepoints
for t1 = 1:length(anaDat)-1
    t2 = t1 + 1;
    qMatrix1 = anaDat(t1).stats.qMatrix;
    qMatrix2 = anaDat(t2).stats.qMatrix;
    
%     if t1 == 18
%         keyboard
%     end
    
    %select Q-matrix for positional uncertainity of tags in tp1, tp2
    Q1=qMatrix1([(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3],...
        [(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3]);
    Q2=qMatrix2([(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3 ],...
        [(tags(1)-1)*3+1:(tags(1)-1)*3+3 (tags(2)-1)*3+1:(tags(2)-1)*3+3]);
    
    Q=blkdiag(Q1,Q2);
    
    %calculate distance vector
    dif1=distanceVector(t1,:);
    dif2=distanceVector(t2,:);
    
    %calculate distance vectors in pixel
    pix2mu = [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z];
    mu2pix = pix2mu.^-1;
    dif1 = dif1.*mu2pix;
    dif2 = dif2.*mu2pix;
    
    %calculate distances
    [dist,normVec] = normList([dif1;dif2]);
    distance1=dist(1,:);
    distance2=dist(2,:);
    sig2muFact(t1,1) = norm(mean(normVec,1).*pix2mu);
    
    %debug
    dis1(t1,1) = distance1;
    dis2(t2,1) = distance2;
    time1(t1) = t1;
    time2(t1) = t2;
    
    %calculate Hessian
    H=[1/distance1*[dif1, -dif1] , -1/distance2*[dif2, -dif2]];
    Qdd=H*Q*H';
    
    %sigmaDistance = noise(at all tag positions) * Qdd
    sigD(t1,1)=sqrt(mean([anaDat(t1).stats.noise(tags(1)),anaDat(t1).stats.noise(tags(2)),...
            anaDat(t2).stats.noise(tags(1)),anaDat(t2).stats.noise(tags(2))])*Qdd);
    
    %test: is distance change significantly larger than sigmaDistance?
    %as the variances for the measured values are a priori unknown, we should
    %use a t-test. The number of degrees of freedom, however, is very large, as
    %all pixels of the image count -> use gauss
    te(t1,1)=abs((distance2-distance1))/sigD(t1); %te has to be positive
    testBound = norminv((1-TEST_PROB/2),0,1);
    if te(t1)>testBound; 
        tRes(t1)=1;
    end;
    
    
end %for-loop

%for error bars: calculate sigmaD in microns
sigD = sigD.*sig2muFact;
