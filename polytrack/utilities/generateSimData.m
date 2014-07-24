nrOfCells = 50;
nrOfFrames = 100;
xSize = 1340;
ySize = 1000;
xCoord = zeros(nrOfCells,nrOfFrames);
yCoord = zeros(nrOfCells,nrOfFrames);
MPM = zeros(nrOfCells,2*nrOfFrames);

% Randomly generate the first set of coordinates
for i = 1 : nrOfCells
    xCoord(i,1) = round(rand(1)*xSize);
    yCoord(i,1) = round(rand(1)*ySize);
end

% Continue by giving each cell a small displacement
for i = 1 : nrOfCells
    for j = 2 : nrOfFrames
        xCoord(i,j) = xCoord(i,j-1) + 1;
        yCoord(i,j) = yCoord(i,j-1) - 1;
        
        % Check for image borders
        if xCoord(i,j) < 1 | xCoord(i,j) > xSize | ...
           yCoord(i,j) < 1 | yCoord(i,j) > ySize     
            xCoord(i,j) = 0;
            yCoord(i,j) = 0;
        end
    end
end

% Store the generated values in an MPM
for i = 1 : nrOfFrames
    MPM(:,2*i-1) = yCoord(:,i);
    MPM(:,2*i) = xCoord(:,i);
end