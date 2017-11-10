
%outDir = uigetdir(pwd,'Select a directory to save the output to:');
outDir = '/home/he19/files/LCCB/nih/test cases/sphere radius sweep';


%Make sure we get all the simple neighborhoods, and more thoroughly sample
%the smaller radii
radii = [0 1 sqrt(2) sqrt(3) linspace(2,10,32) 12:2:24]; 
nRad = numel(radii);

m = cell(nRad,1);

currOutDir = [outDir filesep 'unprocessed noiseless sphere sweep'];
mkClrDir(currOutDir);


padSize = 5;

for j = 1:nRad
            
    %Make a mask with a sphere of this radius in it
    m{j} = binarySphere(radii(j));
    m{j} = padarray(m{j},padSize * ones(1,3));

    try
        %Run the geom analysis
        mProp{j} = analyze3DMaskGeometry(m{j});   
    catch em
        disp(['err on radius ' num2str(radii(j))'])
        mProp{j} = em;
    end
   disp(j/nRad)
       
end

save([currOutDir filesep 'results.mat'],'mProp','m','radii');

jkl=1;