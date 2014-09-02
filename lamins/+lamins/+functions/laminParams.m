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
        case 'ag_080712wt_Reconstructed 3.mat'
            %p.channelOrder = [1 2 3 4];
            p.channels.DAPI = 1;
            p.channels.LaminA = 2;
            p.channels.LaminB1 = 3;
            p.channels.LaminB1B2 = 4;
            p.goodZ = 8:20;
            p.movieNum = 2;
    end
    p.channels.order = [p.channels.DAPI p.channels.LaminA p.channels.LaminB1 p.channels.LaminB1B2];
end
