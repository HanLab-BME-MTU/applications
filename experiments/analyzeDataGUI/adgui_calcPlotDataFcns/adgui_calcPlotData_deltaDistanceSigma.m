function [sigma] = adgui_calcPlotData_deltaDistanceSigma(anaDat,tags,distanceVectors,distances)
%calculate the uncertainity of a distance change (in microns)

%init vars & consts
sigD = zeros(length(anaDat)-1,1);

%loop through all timepoints
for t1 = 1:length(anaDat)-1
    t2 = t1 + 1;
    qMatrix1 = anaDat(t1).stats.qMatrix;
    qMatrix2 = anaDat(t2).stats.qMatrix;
    
    %select Q-matrix for positional uncertainity of tags in tp1, tp2
    Q1=qMatrix1([(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3],...
        [(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3]);
    Q2=qMatrix2([(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3 ],...
        [(tags(1)-1)*3+1:(tags(1)-1)*3+3 (tags(2)-1)*3+1:(tags(2)-1)*3+3]);
    
    Q=blkdiag(Q1,Q2);
    
    %get distance vectors and distances
    distanceV1 = distanceVectors(t1,:);
    distanceV2 = distanceVectors(t2,:);
    distanceN1 = distances(t1);
    distanceN2 = distances(t2);
    
    %calculate Hessian
    H=[1/distanceN1*[distanceV1, -distanceV1] , -1/distanceN2*[distanceV2, -distanceV2]];
    Qdd=H*Q*H';
    
    %sigmaDistance = noise(mean over all tag positions) * Qdd
    sigma(t1,1)=sqrt(mean([anaDat(t1).stats.noise(tags(1)),anaDat(t1).stats.noise(tags(2)),...
            anaDat(t2).stats.noise(tags(1)),anaDat(t2).stats.noise(tags(2))])*Qdd);
    
%     %test: is distance change significantly larger than sigmaDistance?
%     %as the variances for the measured values are a priori unknown, we should
%     %use a t-test. The number of degrees of freedom, however, is very large, as
%     %all pixels of the image count -> use gauss
%     te(t1,1)=abs((distance2-distance1))/sigD(t1); %te has to be positive
%     testBound = norminv((1-TEST_PROB/2),0,1);
%     if te(t1)>testBound; 
%         tRes(t1)=1;
%     end;
    
    
end %for-loop

