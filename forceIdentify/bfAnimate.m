%This script file shows the animation of identified body force over time.

%Load the identified body force.
load([resultPath 'bfId']);

%Get the image of the cell.
rawI = cell(numTimeSteps,1);
for jj = 1:numTimeSteps
   rawI{jj} = imread(imgFile{jj});
end

%M = vectorFieldAnimate([bfDisplayPx bfDisplayPy],recBF,100000, ...
%   'bgImg',rawI{1},'vColor','r','colorMap','default');
M = vectorFieldAnimate([bfDisplayPx bfDisplayPy],recBF,bfScale, ...
   'bgImg',rawI{1},'colorMap','autumn');
