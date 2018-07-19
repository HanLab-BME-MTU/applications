function addRandomNoise(obj,nNoisePoints)
pointsNoise = repmat(obj.data.roiPosition,nNoisePoints,1)+repmat(obj.data.roiSize,nNoisePoints,1).*rand(nNoisePoints,3);
obj.data.points = [obj.data.points;pointsNoise];

disp('Simulation: Random points added!');
end