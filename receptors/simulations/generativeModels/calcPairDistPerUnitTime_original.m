function distValues = calcPairDistPerUnitTime_original(ptcle1,ptcle2,numDispPerUnitTime)
%CALCPAIRDISTPERUNITTIME calculates the distance between two particles
%at a set of points along their path from an initial to a final position.  
%NOTE: this function has been updated (calcPairDistPerUnitTime).
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
%   Robel Yirdaw, 07/07/2014
%
    
    if (ptcle1.xf ~= ptcle1.xi)        
        %Determine steps taken along x
        xStep1 = ptcle1.xf - ptcle1.xi;
        %Calculate incremental displacement sizes along x         
        delX1 = xStep1/numDispPerUnitTime;
        %Construct array of displacements
        xDisp1 = ptcle1.xi:delX1:ptcle1.xf;
    else
        %Incremental displacement along x is constant
        xDisp1(1:(uint32(numDispPerUnitTime)+1)) = ptcle1.xi;
    end
    
    if (ptcle2.xf ~= ptcle2.xi)
        %Determine steps taken along x
        xStep2 = ptcle2.xf - ptcle2.xi;
        %Calculate incremental displacement sizes along x
        delX2 = xStep2/numDispPerUnitTime;
        %Construct array of displacements
        xDisp2 = ptcle2.xi:delX2:ptcle2.xf;
    else
        %Incremental displacement along x is constant
        xDisp2(1:(uint32(numDispPerUnitTime)+1)) = ptcle2.xi;
    end

    if (ptcle1.yf ~= ptcle1.yi)
        %Determine steps taken along y
        yStep1 = ptcle1.yf - ptcle1.yi;
        %Calculate incremental displacement sizes along y
        delY1 = yStep1/numDispPerUnitTime;
        %Construct array of displacements
        yDisp1 = ptcle1.yi:delY1:ptcle1.yf;
    else
        %Construct array of displacements as a constant
        yDisp1(1:(uint32(numDispPerUnitTime)+1)) = ptcle1.yi;
    end
        
    if (ptcle2.yf ~= ptcle2.yi)
        %Determine steps taken along x
        yStep2 = ptcle2.yf - ptcle2.yi;
        %Calculate incremental displacement sizes along y
        delY2 = yStep2/numDispPerUnitTime;    
        %Construct array of displacements
        yDisp2 = ptcle2.yi:delY2:ptcle2.yf;    
    else
        %Construct array of displacements as a constant
        yDisp2(1:(uint32(numDispPerUnitTime)+1)) = ptcle2.yi;
    end
    
    
    %Calculate distances at each displacement point    
    distValues = sqrt( (xDisp2 - xDisp1).^2 + (yDisp2 - yDisp1).^2 );
    
end

