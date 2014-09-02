function steerable = expandSteerable(compressed)
    steerable.res = compressed.res;
    steerable.nms = compressed.res;
    steerable.nms(compressed.nms == 0) = 0;
    steerable.theta = steerable.nms;
    steerable.theta(steerable.nms) = compressed.theta;
end
