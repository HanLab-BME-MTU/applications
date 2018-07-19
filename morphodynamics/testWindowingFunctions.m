function testResults = testWindowingFunctions(testMasks)

%Function for debugging and testing the windowing functions
%Input testMasks should be a cell-array containing logical matrices

showPlots = true;

nMasks = numel(testMasks);

testPerpSizes = [1 10 100];
testParaSizes = testPerpSizes;
nPerpSizes = numel(testPerpSizes);
nParaSizes = numel(testParaSizes);


for j = 1:nMasks        
    for k = 1:nPerpSizes
        for l = 1:nParaSizes
        
            currString = ['Mask ' num2str(j) ',' num2str(testPerpSizes(k)) 'x' num2str(testParaSizes(k)) ':'];
            
            try
                tic;
                testResults(j,k,l).winMat = getMaskWindowsPixelLevel(testMasks{j},testPerpSizes(k),testParaSizes(l));
                testResults(j,k,l).windowingTime = toc;
            catch err
                testResults(j,k,l).error = err;
                disp([currString ' error! = ' err.message]);                                
            end
            
            testResults(j,k,l).perpSize = testPerpSizes(k);
            testResults(j,k,l).paraSize = testParaSizes(l);            
    
        end
    end
end











