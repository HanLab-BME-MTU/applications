function distValues = calcPairDistPerUnitTime(ptcle1,ptcle2,numDispPerUnitTime)
%CALCPAIRDISTPERUNITTIME calculates the distance between two particles
%at a set of points along their path from an initial to a final position.
%
%   Input:
%       1) ptcle1: a struct containing the initial and final coordinates of
%                  particle 1.
%       2) ptcle2: a struct containing the initial and final coordinates of
%                  particle 2.
%       3) numDispPerUnitTime: number of displacement points between the
%                              initial and final positions.
%
%   Output:
%       distValues: the calculated distances between particle 1 and 2 at
%                   each displacement point along their paths.
%
%   Robel Yirdaw, 07/07/14
%       Modified, 10/17/14
%
    
    %Initialize array of displacements as constants
    xDisp1(1:uint64(numDispPerUnitTime)+1) = ptcle1.xi;
    yDisp1(1:uint64(numDispPerUnitTime)+1) = ptcle1.yi;
    xDisp2(1:uint64(numDispPerUnitTime)+1) = ptcle2.xi;
    yDisp2(1:uint64(numDispPerUnitTime)+1) = ptcle2.yi;

    %Determine steps taken along x and y
    xStep1 = ptcle1.xf - ptcle1.xi;
    yStep1 = ptcle1.yf - ptcle1.yi;    
    xStep2 = ptcle2.xf - ptcle2.xi;    
    yStep2 = ptcle2.yf - ptcle2.yi;     

    %Calculate incremental displacement sizes along x and y     
    delX1 = xStep1/numDispPerUnitTime;
    delX2 = xStep2/numDispPerUnitTime;    
    delY1 = yStep1/numDispPerUnitTime;
    delY2 = yStep2/numDispPerUnitTime;    
    
    %Calculate incremental displacement vectors along x and y
    delVecX1 = delX1*(1:numDispPerUnitTime);
    delVecY1 = delY1*(1:numDispPerUnitTime);
    delVecX2 = delX2*(1:numDispPerUnitTime);
    delVecY2 = delY2*(1:numDispPerUnitTime);
    
    %Calculate incremental displacements along x and y
    xDisp1(2:end) = xDisp1(2:end) + delVecX1;
    yDisp1(2:end) = yDisp1(2:end) + delVecY1;
    xDisp2(2:end) = xDisp2(2:end) + delVecX2;
    yDisp2(2:end) = yDisp2(2:end) + delVecY2;       
    
    %Calculate distances at each displacement point    
    distValues = sqrt( (xDisp2 - xDisp1).^2 + (yDisp2 - yDisp1).^2 );
    
end

