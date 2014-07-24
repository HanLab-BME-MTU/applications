function open_all_images

i_am_here = pwd
allpaths = genpath(i_am_here);
marks = [0,strfind(allpaths,':'), length(allpaths)+1];
for t = 1:size(marks,2)-1
    
    localims = dir([allpaths((marks(t)+1):(marks(t+1)-1)),'/*.tif']);
    for i =1:size(localims)
        [allpaths((marks(t)+1):(marks(t+1)-1)), '/',localims(i).name]
        copyfile([allpaths((marks(t)+1):(marks(t+1)-1)), '/',localims(i).name],[i_am_here,'/',num2str(t),num2str(i),'.tif']);
    end
end
    

    