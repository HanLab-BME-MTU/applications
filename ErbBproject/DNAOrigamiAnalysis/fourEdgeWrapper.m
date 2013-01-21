%Script for applying 4Edge filter
%
%

%load(file);
list = extractF(features);
drift = getDrift(features,'amp',2000);
fcorr = correctDrift(features,drift);
list2 = extractF(fcorr);

%removes points used for drift correction
list2(list2(:,3)>2000,:)=[];

n = numel(list2(:,1));

%representive sigma in nanometers
sig = 10.0*2;

%nanometers per pixel
npp = 128.0;

%scaling factor
sf = 10.0;

sigma = (sig/npp)*ones([n,1]);

%makes an image
%img = PointP_Plot2Color([list2(:,[1,2]),sigma,sigma], [0,0,sig,sig], sf, 50,'test.tif');
%[counts,bins]=histND(list2(:,1:2),256*sf);

%Generates bins edges for generating a image by histograming data
inc = 1/sf;
e = 0:inc:256;
edge = [{e},{e}];
[counts,C]=hist3(list2(:,1:2),'Edges',edge);

img = zeros([size(counts),3]);
img(:,:,2)=counts;

%makes and applies masks
m = fourEdgeMask(npp/sf,sig/npp);
f = filterMany(img(:,:,2),m);
f2=f;

n = sum(sum(f > 0));

[d,idx]=sort(f(:),'descend');

% retain only the top 1% of values
f2(idx(ceil(n*0.01):end)) = 0.0;

img(:,:,1)=f2;

imshow(img);

CC = bwconncomp(f2>0);

%The following creates a 2D histogram with 2x scaling of the points in each
%identified object

ObjectRecon = cell([CC.NumObjects,1]);
Area = zeros([CC.NumObjects,1]);
inc2 = 1/(2*sf);
tbd = Area;

for k=1:CC.NumObjects
    %finds the rectange that incloses the object plus a buffer around it
    [i,j] = ind2sub(size(f2),CC.PixelIdxList{k});
    ObjectRecon{k} = struct('counts',[],'bins',[],'PntListIdx',[]);
    %sets a minimum size for reconstruction
    if numel(i) < 3 
        tbd(k)=1;
    end
    
    xmin = edge{1}(min(i))-0.7;
    xmax = edge{1}(max(i))+0.7;
    ymin = edge{2}(min(j))-0.7;
    ymax = edge{2}(max(j))+0.7;
    
    Area(k) = (xmax-xmin-1.4)*(ymax-ymin-1.4); % in pixels
    
    %define bin edges for image generation
    edge2 = [{xmin:inc2:xmax},{ymin:inc2:ymax}];
    
    inWindow = list2(:,1) >= xmin & list2(:,1) <= xmax & list2(:,2) >= ymin & list2(:,2) <= ymax;
    Pnts = list2(inWindow,1:2);
    %[cnt,bn]=hist3(Pnts,[ ceil((xmax-xmin)*sf*2),ceil((ymax-ymin)*sf*2)]);
    [cnt,bn]=hist3(Pnts,edge2);
    ObjectRecon{k}.counts = cnt;
    ObjectRecon{k}.bins = edge2;
    ObjectRecon{k}.PntListIdx = inWindow;
end

CC.ObjectRecon=ObjectRecon;
CC.Area = Area;
ObjectRecon(logical(tbd))=[];
h = ObjectReconViewer(ObjectRecon);