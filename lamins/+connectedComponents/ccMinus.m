function cc = ccPlus(A,B)
    cc = connectedComponents.ccBinaryOp(@setdiff,A,B);
end
