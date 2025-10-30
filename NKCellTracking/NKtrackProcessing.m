%% loading data
clear;
eucDist=@(x,y) sqrt(x^2+y^2);

if usejava('desktop')
    [fileSFolders, pathSFolders] = uigetfile('*.mat','Select MovieData file');
    else
        disp({'Enter Path to MovieData (.mat)';
            ['Your current path: ' pwd]});
        rawPath = input(': ','s');
        if isempty(rawPath)
            pathSFolders = 0;
        else
            [pathSFolders, fileSFolders] = fileparts(rawPath);
        end
end


try 
    MD=MovieData.load(append(pathSFolders,fileSFolders));
catch
    disp(['Error :: failed to load file '  fileSFolders])
end
path=append(pathSFolders,fileSFolders);

iAnalysis=MD.getProcessIndex('MotionAnalysisProcess');
if isempty(iAnalysis)
    error("Motion Anlysis must be run");
end

% motion classes 0:immobile 1:confined Brownian 2:pure Brownian 3:directed motion

analysisProcess=MD.getProcess(iAnalysis);
iAPC=find(analysisProcess.checkChannelOutput);
fp=analysisProcess.outFilePaths_{iAPC};
tr=load(fp); %tr contains 'tracks' & 'diffAnalysisRes'
X=MD.reader.getSizeX;
Y=MD.reader.getSizeY;
Z=MD.reader.getSizeZ;
T=MD.reader.getSizeT;
ps=MD.pixelSize_; %pixels size in nanometres 
if Z>1
    error("Image must be 2D");
    %implement stacking
end
img=zeros(X,Y);
vImg=img;
cImg=img;
for i=1:T
    im=MD.reader.loadImage(1,i);
    fImg=vesselFilter2(im);
    vImg=logical(vImg+fImg);
    im=MD.reader.loadImage(2,i);
    fImg=cellFilter(im);
    cImg=logical(cImg+fImg);

end
disp(['Finished Creating Image Masks']);
%%  merge data
tracks=tr.tracks;
[numTracks,~]=size(tracks);
rel = repmat(struct('x',[],'y',[],'in',[],'c',[],'distance',[]),numTracks,1);
cl = repmat(struct('class',[],'mssSlope',[],'diffCoef',[]),numTracks,1);
locVel=cell(3,1); %in vessel, out cancer out vessel, in cancer & out vessel
locMSS=cell(3,1);
locCoef=cell(3,1);
class={tr.diffAnalysisRes.classification};
dim={tr.diffAnalysisRes.fullDim};
clCount=[0,0,0,0];


for i = 1:numTracks
    dT=size(tracks(i).tracksFeatIndxCG,2);
    t=tracks(i).tracksCoordAmpCG;
    [n,~]=size(t);
    x=nan(n,dT);
    y=nan(n,dT);
    b=false(n,dT);
    c=false(n,dT);
    dis=zeros(n,dT);
    classification=nan(n,1);
    mss=nan(n,1);
    diffCoef=nan(n,1);



    for j=1:n


        classification(j)=class{i}(j,2);
        if ~isnan(classification(j))
            clCount(classification(j)+1)=clCount(classification(j)+1)+1;
        end

        rT=reshape(t(j,:),8,dT);
        x(j,:)=rT(1,:);
        y(j,:)=rT(1,:);
        for k=1:dT
            lX=round(x(j,k));
            lY=round(y(j,k));
            if isnan(lX) || isnan(lY) || isnan(classification(j))
                continue
            else 
                b(j,k)=vImg(lX,lY);
                c(j,k)=cImg(lX,lY);
                if k~=1 %&& classification(j)
                    dis(j,k)=eucDist((x(j,k)-x(j,k-1)),(y(j,k)-y(j,k-1)))*ps; %nanometres per interval
                    if b(j,k) %if in vessel
                        locVel{1}=[locVel{1},dis(j,k)];
                    else 
                        if c(j,k) %if out of vessel in cancer
                            locVel{3}=[locVel{3},dis(j,k)];
                        else    %if out of cessel out of cancer
                            locVel{2}=[locVel{2},dis(j,k)];
                            
                        end

                    end
                    
                end
            end
        end
        

        mss(j)=dim{i}.mssSlope(j);
        diffCoef(j)=dim{i}.normDiffCoef(j);
        mb=mode(b(j,:));
        mc=mode(c(j,:));
        if mb %if in vessel
            locMSS{1}=[locMSS{1},mss(j)];
            locCoef{1}=[locCoef{1},diffCoef(j)];
        else
            if mc %if in cell
                locMSS{3}=[locMSS{3},mss(j)];
                locCoef{3}=[locCoef{3},diffCoef(j)];

            else %if not in cell
                locCoef{2}=[locCoef{2},diffCoef(j)];
                locMSS{2}=[locMSS{2},mss(j)];

            end
            
        end
    end

 

    rel(i).x=x;%rT(1,:);
    rel(i).y=y;%rT(2,:);
    rel(i).in=b;
    rel(i).distance=dis;
    rel(i).c=c;
    cl(i).class=classification;
    cl(i).mssSlope=mss;
    cl(i).diffCoef=diffCoef;


end

figure();
boxPlotCellArray(locVel',{'V','~V ~C', '~V C'});
title("Velocity");
figure();
boxPlotCellArray(locMSS',{'V','~V ~C', '~V C'});
title("MSS");
figure();
boxPlotCellArray(locCoef',{'V','~V ~C', '~V C'});
title("coef");
figure();
bLabel=categorical({'Immobile','Confined Brownian','Pure Brownian','Directed Motion'});
bLabel=reordercats(bLabel,{'Immobile','Confined Brownian','Pure Brownian','Directed Motion'});
bar(bLabel,clCount);
title("Movement Classification")
%img=MD.reader.loadImage(1,1)





