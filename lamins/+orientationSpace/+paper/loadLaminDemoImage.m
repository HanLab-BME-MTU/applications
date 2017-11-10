function [ I, MD ] = loadLaminDemoImage(  )
%loadLaminDemoImage Load demo lamin image on multiple platforms

if(~exist('I','var'))
    [~,hostname] = system('hostname');
    hostname = strtrim(hostname);

    switch(hostname)
        case 'FSM2UA220003Q'
            % HP Z820 Goldman workstation
            cd 'Z:\Takeshi\N-SIM\040715';
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.nd2');
            I = MD.channels_(1).loadImage(1,10);
        case 'mkitti-jaqaman'
            % Laptop T440s
            load('C:\Users\Mark Kittisopikul\Documents\Data\Lamins\MEFLB1-LACLB12-006_Reconstructed_study\MEFLB1-LACLB12-006_Reconstructed\MEFLB1-LACLB12-006_Reconstructed.mat');
            MD.sanityCheck;
            I = MD.channels_(1).loadImage(1,10);
        case 'FSMPC0KTM9U'
            cd 'P:\Basic_Sciences\CMB\GoldmanLab\Takeshi\N-SIM\040715';
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.nd2');
            I = MD.channels_(1).loadImage(1,10);
        otherwise
            % BioHPC
            cd ~/shortcuts/MEFLB1-LACLB12-006_Reconstructed/
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.mat');
            I = MD.channels_(1).loadImage(1,10);
    end
end



end

