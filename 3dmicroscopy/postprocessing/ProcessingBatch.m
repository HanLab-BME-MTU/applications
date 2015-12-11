%%

%script to batch process binary data

clc; clear all; close all

folder='//project/cellbiology/gdanuser/microscopeDevelopment/multiplane/RPE/tractinGFP/151208/Cell7/'


Xdim=512;
Ydim=128
Zdim=68

%import binary files and convert to 3DStack

im1=BinaryTo3D([folder '1_CAM01_000001.bin'],Xdim,Ydim,Zdim);
im2=BinaryTo3D([folder '1_CAM02_000001.bin'],Xdim,Ydim,Zdim);
im3=BinaryTo3D([folder '1_CAM03_000001.bin'],Xdim,Ydim,Zdim);

%% maximum intensity projection
pro1=squeeze(max(permute(im1,[2,3,1])));
pro2=squeeze(max(permute(im2,[2,3,1])));
pro3=squeeze(max(permute(im3,[2,3,1])));

%figure;imshow(pro1,[]); figure;imshow(pro2,[]);figure;imshow(pro3,[])

%% stitch images

stitchedView=StitchProjections(pro1,pro2,pro3,[31,-7;-7,11]);

figure;imshow(stitchedView,[])