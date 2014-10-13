function p = laminParams(MD)
    p.channels.labels = {'DAPI', 'Lamin A','Lamin B1','anti-Lamin B1, B2'};
    switch(MD.getFilename())
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
    end
    p.channels.order = [p.channels.DAPI p.channels.LaminA p.channels.LaminB1 p.channels.LaminB1B2];
end
