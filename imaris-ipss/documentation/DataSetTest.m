% connect to Imaris
vImarisApplication = actxserver('Imaris.Application');

% load dataset from file
vImarisApplication.FileOpen('Q:\Bitplane\Imaris4\images\wubbel.ims');
vImarisApplication.mTimePoint = 3;

% get and display the middle slice of current data set (channel 0 timepoint 0)
vImarisDataSet = vImarisApplication.mDataSet;
vSizeZ = vImarisDataSet.mSizeZ;
vDataSetMiddleSlice = vImarisDataSet.GetDataSlice(vSizeZ / 2, 1, 0);
figure,imshow(vDataSetMiddleSlice);

% make a new data set
vDataSet = uint8(255*rand(256,256,32,1,4));
vImarisDataSet.SetData(vDataSet);
