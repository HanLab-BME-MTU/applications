% make fake movies

% rand velocity (1)
% rand size -/+ 5 (1)
% rand direction (every iteration)

nFrames = 20;
nMovies = 150;
imSize = [200 200];
noiseMean = .01;
noiseVar = .005;

movies = cell(nMovies, 1);

for i = 1:nMovies
    
    obj = genCell;
    imStack = zeros(imSize(1,1), imSize(2), nFrames);
    
    for fi = 1:nFrames
        I = GaussMask2D(obj.sizeCell, imSize, [obj.ypos, obj.xpos]);
        I = imnoise(I,'gaussian',noiseMean, noiseVar);
        obj = updatePos(obj);
        imStack(:,:,fi) = I;
    end

    movies(i) = {imStack};
end

disp('done');

% playMovie(movies{1})

function playMovie(movie)

    figure();
    for i=1:size(movie,3)
        imshow(movie(:,:,i));
        input('press enter');
    end

end


function [obj] = genCell

    obj= struct();
    obj.xpos = randi([-50 0]);
    obj.ypos = randi([-50 0]);
    obj.xvel = randi([0 7]);
    obj.sizeCell = randi([20 50]);

end

function [obj] = updatePos(obj)

    obj.xpos = obj.xpos + obj.xvel;
    obj.ypos = obj.ypos + randi([-3 3]);
    obj.sizeCell = obj.sizeCell + randi([-2 2]);
    
end

% xvel will be random