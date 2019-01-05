function out = Integral_mu_inf(mu,chi,xi)
%
tempA1 = exp(-mu*chi)/chi;
tempA2 = xi*igamma(0,mu*chi);
tempA3 = 0;
for uu=2:2
    tempB1 = ((-1)^uu)*(xi^uu)/(factorial(uu));
    tempB21 = exp(-(mu*chi));
    tempB22 = 0;
    for vv = 1:(uu-1)
        temp = (factorial(vv-1))*((-chi)^(uu-vv-1))/...
            ((factorial(uu-1))*(mu^vv));
        tempB22 = tempB22 + temp;
    end
    tempB23 = ((-chi)^(uu-1))/(factorial(uu-1))*(ei(-mu*chi));
    tempB2 = tempB21*tempB22-tempB23;
    tempA3 = tempA3 + tempB1*tempB2;
    
end
out = tempA1-tempA2+tempA3;