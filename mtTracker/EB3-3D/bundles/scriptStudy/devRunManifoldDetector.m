earlyLateML=MovieList.load('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\analysis\earlyAndLateML.mat');
%%
for i=1:1%earlyLateML.getSize()
    MD=earlyLateML.getMovie(i);
    %pack=run1DManifoldDetectorPackage(MD,'package',MD.getPackage(333));
    pack=run1DManifoldDetectorPackage(MD);
    MD.setPackage(333,pack)
    MD.save();
end

