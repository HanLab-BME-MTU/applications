function p = laminParams(MD)
    filename = MD.getFilename();
    % translate equivalent images to original filenames
    switch(filename)
        case 'ag_080712wt_Z15_C3.mat'
            filename = 'ag_080712wt_Reconstructed 3.mat';
    end
    p.channels.labels = {'DAPI', 'Lamin A','Lamin B1','anti-Lamin B1, B2'};
    p.channels.DAPI = 0;
    p.channels.LaminA = 0;
    p.channels.LaminB1 = 0;
    p.channels.LaminB1B2 = 0;
    p.goodZ = 1;
    p.movieNum = -1;
    p.steerable.sigma = 2;

    switch(filename)
        case 'ag_072612_wt_Reconstructed 2.mat'
            %p.channelOrder = [3 1 2 4];
            p.channels.DAPI = 3;
            p.channels.LaminA = 1;
            p.channels.LaminB1 = 2;
            p.channels.LaminB1B2 = 4;
            p.goodZ = 12:28;
            p.movieNum = 1;
            p.steerable.sigma = 5;
        case 'ag_080712wt_Reconstructed 3.mat'
            %p.channelOrder = [1 2 3 4];
            p.channels.DAPI = 1;
            p.channels.LaminA = 2;
            p.channels.LaminB1 = 3;
            p.channels.LaminB1B2 = 4;
            p.goodZ = 8:20;
            p.movieNum = 2;
            p.steerable.sigma = 5;
        case 'Reconstructed TIRF SIM-01.mat'
            p.channels.DAPI = 0;
            p.channels.LaminA = 0;
            p.channels.LaminB1 = 0;
            p.channels.LaminB1B2 = 0;
            p.goodZ = 1;
            p.movieNum = 3;
            p.steerable.sigma = 2;
        case 'Reconstructed TIRF SIM-02.mat'
            p.channels.DAPI = 0;
            p.channels.LaminA = 0;
            p.channels.LaminB1 = 0;
            p.channels.LaminB1B2 = 0;
            p.goodZ = 1;
            p.movieNum = 4;
            p.steerable.sigma = 2;
        case 'MEFLB1M20-01constructed.mat'
            p.movieNum = 5;
        case 'MEFLB1M20-10constructed.mat'
            p.movieNum = 6;
        case 'MEFWTLA-05constructed.mat'
            p.movieNum = 7;
        case 'MEFWTLA-09constructed.mat'
            p.movieNum = 8;
        otherwise
            disp(['Unknown movie parameters: ' filename]);
    end
    p.channels.order = [p.channels.DAPI p.channels.LaminA p.channels.LaminB1 p.channels.LaminB1B2];
end
