function cc = ccPlus(A,B)
    cc = connectedComponents.ccBinaryOp(@union,A,B);
end
