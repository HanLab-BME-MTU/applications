function tRes=testSigSeparationChange(idlist,dataProperties,timePt,tags)
%calculate whether a distance change is significant knowing the uncertainities
%of the positional estimates (0 if not, +1 if distance increases, -1 if distance decreases)

%init vars & consts
tRes=0;
T_TEST_PROB=0.05;

%select Q-matrix for positional uncertainity of tags in tp1, tp2
if ~isempty(idlist(timePt(1)).info.trackQ_Pix)
    Q1=idlist(timePt(1)).info.trackQ_Pix([(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3],...
        [(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3]);
else
    Q1=idlist(timePt(1)).info.detectQ_Pix([(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3 ],...
        [(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3]);
end;
if ~isempty(idlist(timePt(2)).info.trackQ_Pix)
    Q2=idlist(timePt(2)).info.trackQ_Pix([(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3 ],...
        [(tags(1)-1)*3+1:(tags(1)-1)*3+3 (tags(2)-1)*3+1:(tags(2)-1)*3+3]);
else
    Q2=idlist(timePt(2)).info.detectQ_Pix([(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3 ],...
        [(tags(1)-1)*3+1:(tags(1)-1)*3+3, (tags(2)-1)*3+1:(tags(2)-1)*3+3]);
end;
Q=blkdiag(Q1,Q2);

%calculate distance vector
dif1=(idlist(timePt(1)).linklist(tags(1),9:11)-idlist(timePt(1)).linklist(tags(2),9:11));
dif2=(idlist(timePt(2)).linklist(tags(1),9:11)-idlist(timePt(2)).linklist(tags(2),9:11));

%calculate distance vectors in pixel
mu2pix = [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z;...
    dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z].^-1;
dif1 = dif1.*mu2pix;
dif2 = dif2.*mu2pix;

%calculate distances
distance1=sqrt(sum(dif1.^2));
distance2=sqrt(sum(dif2.^2));

%calculate Hessian
H=[1/distance1*[dif1, -dif1] , -1/distance2*[dif2, -dif2]];
Qdd=H*Q*H';

%sigmaDistance = noise * Qdd
sigD=sqrt(mean(idlist(timePt(1)).info.noise)*Qdd);

%test: is distance change significantly larger than sigmaDistance?
te=abs((distance2-distance1))/sigD;
if te>tinv(1-(T_TEST_PROB/2),1);
    tRes=sign(sqrt(sum((dif2.^2)))-sqrt(sum((dif1.^2))));
end;
