%This script file shows the animation of identified body force over time.

%Load the identified body force.
load([resultPath 'bfId']);

%Get the image of the cell.
rawI = cell(numTimeSteps,1);
for jj = 1:numTimeSteps
   rawI{jj} = imread(imgFile{jj});
end

M = vectorFieldAnimate([bfDisplayPx bfDisplayPy],recDispU,dispScale, ...
   'bgImg',rawI{1},'colorMap','autumn');
