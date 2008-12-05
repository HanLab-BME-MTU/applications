function popVelScatter(projData)

close all

for p=1:length(projData)

    homeDir=pwd;

    proj=projData(p,1);

    cd(proj.anDir)

    % get list of files and image size
    [listOfImages] = searchFiles('.tif',[],proj.imDir,0);
    fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
    img = double(imread(fileNameIm));
    [imL,imW] = size(img);

    velLimit=40;
    cMapLength=128;
    cMap=jet(cMapLength);

    temp=proj.trackVelocities;
    m=max(abs([min(temp(:)); max(temp(:))]));
    m=min(m,velLimit);
    temp(temp<-m)=-m;
    temp(temp>m)=m;

    mapper=linspace(-m,m,cMapLength)';

    scrsz = get(0,'ScreenSize');
    screenW=scrsz(3);
    screenL=scrsz(4);
    magCoef=inf;

    maxMagCoefW = (0.8*screenW)/(imW);
    maxMagCoefL = (0.8*screenL)/(imL);
    calcMagCoef = min([magCoef; maxMagCoefW; maxMagCoefL]);

    movieL = round(calcMagCoef*(imL));
    movieW = round(calcMagCoef*(imW));

    figure('Position',[round(screenW*(1-movieW/screenW)/2) round(screenL*(1-movieL/screenL)/2) movieW movieL])

    eval(['MakeQTMovie start ', 'velDistMovie' '.avi']);


    for i=1:20 %size(listOfImages,1)


        fileNameIm = [char(listOfImages(i,2)) filesep char(listOfImages(i,1))];
        img = double(imread(fileNameIm));
        [imL,imW] = size(img);

        imagesc(img)
        colormap gray

        hold on

        D=createDistanceMatrix(temp(:,i),mapper);
        [sD,idx]=sort(abs(D),2);
        scatter(proj.xCoord(:,i),proj.yCoord(:,i),'Marker','.','cData',cMap(idx(:,1),:));

        text(40,40,['Frame: ' num2str(i), ', velLimit = ' num2str(m)],'Color','w','FontWeight','bold','HorizontalAlignment','left')
        MakeQTMovie addaxes
        MakeQTMovie('framerate', 5);
    end

    MakeQTMovie finish

    close all

    cd(homeDir)
end
