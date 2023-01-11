function speedFrequencySave(varargin)
    p=inputParser;
    addOptional(p,'path',0,@ischar)
    parse(p,varargin{:})
    importedPath=p.Results.path;
    
    if ~importedPath    
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
    else
        try 
            MD=MovieData.loadMatFile(importedPath);
        catch
            disp(['Error :: failed to load imported path ' importedPath ])
        end
        path=importedPath;
    end
    disp(['Loaded MovieData ::' MD.movieDataFileName_])
    disp(['Analysing...'])
    
    [speedCell,protSpeedCell] = quantifyMovieFlowSpeed(MD);
    
    %speedCell (window index, depth, time)
    
    
    speed=squeeze(speedCell(:,1,:));
    r=input("Range:");
    speed=speed(r,:);
    [m, n]=size(speed);
    
    num=20;%m;
    fs=1/(6);
    
    %speed=detrend(speed,2);
    
    % figure()
    % title(['Speed ' num2str(num) ' segments'])
    % subtitle([ti]);
    % xlabel('Time (seconds)');
    % hold on
    % for i=[1:20];%floor(linspace(1,m,20))
    % plot([1:length(speed(i,:))]*6,speed(i,:));
    % end
    % hold off
    
    %%FFT
    % figure()
    % title(['Power Spectrum ' num2str(num) ' segments']);
    % subtitle([ti]);
    % xlabel("Frequency (Hz)");
    % ylabel("Power")
    % hold on
    %tiledlayout(floor((length(speed(i,:))+1)/4),4)
    for j=1;%floor(linspace(1,m,floor(m/20))) 
    %     figure()
    %     title(['Power Spectrum ' num2str(j) '-' num2str(j+20) ' segments']);
    %     subtitle([ti]);
    %     xlabel("Frequency (Hz)");
    %     ylabel("Power")
    %     hold on
    %     peaks=cell(20,1);
        for i=[1:length(speed)];%[j:j+20]%floor(linspace(1,m,20))
        %    nexttile
            l=length(speed(i,:));
            ff=fft(speed(i,:)-mean(speed(i,:)));
            p=abs(ff).^2;
            p1=p(1:floor((n+1)/2));
            f=fs/l*(0:floor((l-1)/2));
            %plot(f,p1);
            
            overIndex=find(round(f,4)==round(1/(31*6),4));
            p1(overIndex)=0;
            %maxs=islocalmax(p1);
            %maxs(overIndex)=0;
            [m,maxs]=max(p1);
            maxf=f(maxs);
            peaks{i}=maxf;
            ps{i}={f,p1};
        end
    end
    
    box=[];
    for i=[1:length(speed)]
        for j=[1:length(peaks{i})]
            box=[box,peaks{i}(j)];
            
        end
    end
    speedOut=speed;
    
    name=input("Name:",'s');

    
    
    
    save(string([name,'.mat']),'speedOut','box','path','ps');
end