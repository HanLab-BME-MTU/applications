%% Prelim
path = '/home/sb234/Desktop/CSUX560_Adcon_1a_062609/ch488';
imageFiles = dir([path filesep 'roi' filesep '*.tif']);
h = ones(3);
h(2, 2) = 0;
%% Computation
for i = 1:length(imageFiles)
    I = imread([path filesep 'roi' filesep imageFiles(i).name]);
    Irec = imWTATrouDenoising(I, 3, 0, 3);
    BW = zeros(size(I), 'uint8');
    BW(uint16(Irec) ~= 0) = 1;
    BWfiltered = imfilter(BW, h);
    BW(BWfiltered == 0) = 0;
    imwrite(logical(BW), [path filesep 'analysis' filesep 'masks' ...
        filesep imageFiles(i).name], 'Compression', 'none');
end