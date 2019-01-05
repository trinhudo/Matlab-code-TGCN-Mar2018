function out = TriNhuAppoxIntegral(a,b,c)
%
out = 0;
for nn=0:10
    B1 = ((-1)^nn)*(c^nn)/(factorial(nn))/(b^(-1-nn));
    B2 = igamma(-1-nn,b/a);    
    B = B1*B2;
    out = out + B;
end