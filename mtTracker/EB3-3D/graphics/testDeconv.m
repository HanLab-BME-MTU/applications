addpath(genpath('C:\Users\Philippe\code\utsw-ssh'))
load('C:\Users\Philippe\project-local\dataManagement\testPrintMIPArray\analysis\earlyAndLateML.mat')
MD=ML.getMovies
MD=MD{1}
PSF=stackRead('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\analysis\mbPSF_488_NAp5nap42_z100nm.tif');
vol=MD.getChannel(1).loadStack(1);
XYMIP=(computeMIPs(vol,2.16,0,100));
ZRatio=2.16;

%% 2D deconv with very approximate 2D PSF
testSliceNB=20;
testSliceNBISO=round(testSliceNB*ZRatio);
slice=vol(:,:,testSliceNB);
decslice=deconvlucy(double(slice),double(PSF(:,:,60)),15,0, ones(size(slice)),0,1);

imshow([mat2gray(slice) mat2gray(imresize(decslice,size(slice,1)/size(decslice,1)))])


%% 3D deconv
% resize in Z. 
inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
    [1 MD.getDimensions('X')],[1 MD.getDimensions('Y')],[1 MD.getDimensions('Z')*ZRatio]);
outputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') ceil(MD.getDimensions('Z')*ZRatio)], ...
    [1 MD.getDimensions('X')],[1 MD.getDimensions('Y')],[1 MD.getDimensions('Z')*ZRatio]);

volISO=imwarp(vol,inputRef,affine3d(),'OutputView',outputRef);
%
%volISO = edgetaper(double(volISO),double(PSF));
decvol=deconvlucy(double(volISO),double(PSF),10,0, ones(size(volISO)),0,1);
%decvol=deconvlucy(double(volISO),double(PSF),15);
%estimated_nsr = 100 / var(double(volISO(:)));
%decvol=deconvwnr(double(volISO),double(PSF),estimated_nsr);

imshow([mat2gray(slice) ...
        mat2gray(volISO(:,:,testSliceNBISO)) ...
        mat2gray(imresize(decslice,size(slice,1)/size(decslice,1))) ]);
    imshow( [mat2gray(imresize(decvol(:,:,testSliceNBISO),size(slice,1)/size(decvol(:,:,testSliceNBISO),1))) ])

%% trial on non-deskewed data (senseless shit)
tmp=load('C:\Users\Philippe\project-local\externBetzig\raw\adavid\smallSamples\A1_HeLa_Cells_EB1\prometaphase\earlyCell1_12_raw\analysis\movieDataPixelSize.mat');
MD=tmp.MD.addAnalysisFolder('C:\Users\Philippe\project-local\externBetzig\raw\adavid\smallSamples\A1_HeLa_Cells_EB1\','C:\Users\Philippe\project-local\externBetzig\analysis\testDeconv');
vol=MD.getChannel(1).loadStack(1);
slice=vol(:,:,testSliceNB);
ZRatio=4;

inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
    [1 MD.getDimensions('X')],[1 MD.getDimensions('Y')],[1 MD.getDimensions('Z')*ZRatio]);
outputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') ceil(MD.getDimensions('Z')*ZRatio)], ...
    [1 MD.getDimensions('X')],[1 MD.getDimensions('Y')],[1 MD.getDimensions('Z')*ZRatio]);

volISO=imwarp(vol,inputRef,affine3d(),'OutputView',outputRef);
decvol=deconvlucy(double(volISO),double(PSF),20,0, ones(size(volISO)),0,1);
%decslice=deconvlucy(double(slice),double(PSF(:,:,60)),20,0, ones(size(slice)),0,1);

imshow([mat2gray(slice) ...
        mat2gray(volISO(:,:,round(testSliceNB*ZRatio)))
        mat2gray(imresize(decslice,size(slice,1)/size(decslice,1))) ...
        mat2gray(imresize(decvol(:,:,round(testSliceNB*ZRatio)),size(slice,1)/size(decvol(:,:,round(testSliceNB*ZRatio)),1))) ])