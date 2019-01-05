function [pN,pF] = PowerAllocation(RthN,RthF)
pN = (2^(2*RthN) - 1)/(2^(2*RthN+2*RthF)-1);
pF = 1-pN;
