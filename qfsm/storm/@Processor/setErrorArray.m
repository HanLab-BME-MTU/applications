function setErrorArray(obj,errorX,errorY,errorZ)

obj.data.error = [ones(obj.data.nPoints,1)*errorX ...
    ones(obj.data.nPoints,1)*errorY ...
    ones(obj.data.nPoints,1)*errorZ];

end
