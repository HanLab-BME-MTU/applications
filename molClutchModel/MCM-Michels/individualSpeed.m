clc
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

disp(['Loaded MovieData ::' MD.movieDataFileName_])
disp(['Analysing...'])

ti=strrep(fileSFolders,'_','\_');
[speedCell,protSpeedCell] = quantifyMovieFlowSpeed(MD);

%speedCell (window index, depth, time)


speed=squeeze(speedCell(:,1,:));
[m, n]=size(speed);

num=m;%20;
fs=1/(6);

%speed=detrend(speed,2);

figure()
title(['Speed ' num2str(num) ' segments'])
subtitle([ti]);
xlabel('Time (seconds)');
hold on
for i=floor(linspace(1,m,20))
plot([1:length(speed(i,:))],speed(i,:));
end
hold off

%%FFT
figure()
title(['Power Spectrum ' num2str(num) ' segments']);
subtitle([ti]);
xlabel("Frequency (Hz)");
ylabel("Power")
hold on
%tiledlayout(floor((length(speed(i,:))+1)/4),4)
for j=floor(linspace(1,m,floor(m/20))) 
    figure()
    title(['Power Spectrum ' num2str(j) '-' num2str(j+20) ' segments']);
    subtitle([ti]);
    xlabel("Frequency (Hz)");
    ylabel("Power")
    hold on
    for i=[j:j+20]%floor(linspace(1,m,20))
    %    nexttile
        l=length(speed(i,:));
        ff=fft(speed(i,:));
        p=abs(ff).^2;
        p1=p(1:floor((n+1)/2));
        p1(1)=0;
        f=fs/l*(0:floor((l-1)/2));
        plot(f,p1);
    
    end
end
