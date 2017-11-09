function processCell = queryMovieListProcess(ML,name,varargin)
% load scoring daip = inputParser;
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('select','all');
ip.parse(varargin{:});
p=ip.Results;
processCell=cell(ML.getSize(),1);
movieCell=ML.getMovies();

for MDIdx=1:length(movieCell)
    MD=movieCell{MDIdx};
    proc=[];
    try
    proc=MD.findProcessTag(name,'queryField','name');
    catch 
    end;
    switch p.select
        case 'all'
            processCell{MDIdx}=proc; 
        case 'last'
            processCell{MDIdx}=proc(end);
    end
end