function sigma = adgui_calcPlotData_distanceSigma(anaDat,tags,distanceVectors,distances,centers,deltaT);
%calculates the uncertainity of a distance measurement using gaussian error
%propagation in adgui_calcPlotData (all in microns)

%init vars & consts
sigD = zeros(size(distanceVectors,1),1);

distances(find(distances == 0))=NaN;

%check if optional deltaT has been supplied
if nargin<6|isempty(deltaT)
    deltaT = [];
end

switch length(tags)+10*~isempty(deltaT)
    
    case 2 %we calculate a distance
        
    deltaTime = 0;
    Q3 = zeros(3);
    Q4 = zeros(3);
    doQ = 0;
    
case 1 %we calculate a displacement
    
    deltaTime = 1;
    tags(2) = tags(1);
    
    if ~centers %uncorrected displacement
        Q3 = zeros(3);
        Q4 = zeros(3);
        doQ = 0;
    else
        doQ = 1;
    end
    
    case 11 %diffusion
    
    deltaTime = deltaT;
    tags(2) = tags(1);
    
    if ~centers %uncorrected displacement
        Q3 = zeros(3);
        Q4 = zeros(3);
        doQ = 0;
    else
        doQ = 1;
    end
        
otherwise
    
    error('bad number of tags')
    
end

%loop through all timepoints
for t1 = 1:length(anaDat)-deltaTime
    t2 = t1+deltaTime;
    qMatrix1 = anaDat(t1).stats.qMatrix;
    qMatrix2 = anaDat(t2).stats.qMatrix;
    
    %select Q-matrix for positional uncertainity of the tags
    Q1=qMatrix1([(tags(1)-1)*3+1:(tags(1)-1)*3+3],...
        [(tags(1)-1)*3+1:(tags(1)-1)*3+3]);
    Q2=qMatrix2([(tags(2)-1)*3+1:(tags(2)-1)*3+3 ],...
        [(tags(2)-1)*3+1:(tags(2)-1)*3+3]);
    
    if doQ
        Q3=qMatrix1([(centers(1)-1)*3+1:(centers(1)-1)*3+3],...
            [(centers(1)-1)*3+1:(centers(1)-1)*3+3]);
        Q4=qMatrix2([(centers(1)-1)*3+1:(centers(1)-1)*3+3 ],...
            [(centers(1)-1)*3+1:(centers(1)-1)*3+3]); 
    end
    
    Q=blkdiag(Q1,Q2,Q3,Q4);
    
    %get right distance vector and the corresponding norm (in microns)
    distanceV = distanceVectors(t1,:);
    distanceN = distances(t1);
 
    
    %calculate Hessian in mu
    H = [1/distanceN*[-distanceV, distanceV, distanceV, -distanceV] ];
    Qdd = H*Q*H';
    
    %sigmaDistance = noise * Qdd
    sigma(t1,1) = sqrt(mean([anaDat(t1).stats.noise(tags(1)),anaDat(t2).stats.noise(tags(2))])*Qdd);
    
end %for-loop