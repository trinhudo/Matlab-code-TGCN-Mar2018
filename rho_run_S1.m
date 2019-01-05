close all
clear all
%
%% Simulation parameters
%
K       = 3;                      % # of antenna
rho     = .1:.01:.9;         % power splitting ratio
alpha   = 0;         % time fraction for EH
PS_dB   = 0;                % transmit SNR = Ps/N0 in dB
PS      = 10.^(PS_dB./10);
naN     = (10^(-7))*1e6;    % naN = -100 dBm, BW = 1 MHz
ncN     = (10^(-6))*1e6;    % naN = -90 dBm,  BW = 1 MHz
naF     = (10^(-7))*1e6;
ncF     = (10^(-6))*1e6;
epsilon = 3;                % pathloss exponent
dSF     = 10;                   % S-F distance
dSN     = 3;
dNF     = dSF - dSN;
L       = 1e3;                  % path-loss at reference distance
%
lSN     = L*dSN^-3;             % lambda
lSF     = L*dSF^-3;
lNF     = L*dNF^-3;
%
eta     = 0.7;              % energy conversion coefficient
RthN    = .1;                % target data rate of User N bits/s/Hz
RthF    = .1;               % target data rate of User N bits/s/Hz
[pN,pF] = PowerAllocation(RthN,RthF);
%
SimTimes = 10^0;            % Monte-Carlo repetitions
%
%% Simulation
%
for ss = 1:length(PS_dB)
    for aa = 1:length(alpha)
        fprintf('alpha = %d \n',alpha(aa))
        for rr = 1:length(rho)
            fprintf('rho = %d \n',rho(rr))
            %
            g2 = 2^(RthF*2/(1-alpha(aa))) - 1; % gamma_2 for User F
            % Channel modelling
            for ii = 1:K
                hiF(:,ii) = sqrt(lSF/2)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
                hiN(:,ii) = sqrt(lSN/2)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            end
            hNF = sqrt(lNF/2)*...
                (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            % Channel gains
            giN     = abs(hiN.^2);
            giF     = abs(hiF.^2);
            gNF     = abs(hNF.^2);
            % SNRs
            snriNxF = (1-rho(rr)).*pF.*PS(ss).*giN./...
                ((1-rho(rr)).*pN.*PS(ss).*giN ...
                + (1-rho(rr))*naN + ncN);
            snriNxN = (1-rho(rr)).*pN.*PS(ss).*giN/...
                ((1-rho(rr))*naN + ncN);
            snriF   = pF.*PS(ss).*giF./(pN.*PS(ss).*giF + naF + ncF);
            snrNF   = eta.*PS(ss).*giN.*gNF.*...
                (2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
            %
            % Find the best antenna for User F based on end-to-end SNR
            snrFe2e_i(:,:) = min(snriNxF(:,:),max(snrNF(:,:),snriF(:,:)));            %             end
            [snrFe2e_b,I] = max(snrFe2e_i,[],2);
            % count outage events
            % method 2
            count_1 = snrFe2e_b < g2;
            OP_S1_F_sim(aa,rr) = sum(count_1)/SimTimes;
            %% Analysis
            a1 = (1-rho(rr))*pF*PS(ss)/((1-rho(rr))*naN + ncN);
            a2 = (1-rho(rr))*pN*PS(ss)/((1-rho(rr))*naN + ncN);
            b1 = pF * PS(ss) / (naF + ncF);
            b2 = pN * PS(ss) / (naF + ncF);
            c  = eta*PS(ss)*(2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
            mu_a = g2/(a1-a2*g2);
            mu_b = g2/(b1-b2*g2);
            %
            term1 = 0;
            for ii=0:1:K
                for jj=(K-ii):-1:0
                    kk = K - (ii+jj);
                    A1 = factorial(K)/factorial(ii)/factorial(jj)/factorial(kk);
                    A2 = (1-exp(-mu_a/lSN-mu_b/lSF))^ii;
                    A3 = ((-1)^jj)*exp(-kk*mu_b/lSF)/lNF;
%                     A5 = integral(fun1,0,inf);
                    chi = 1/lNF;
                    if (jj+kk)==0
                        A5 = lNF;
                    else
                    A5 = sqrt(4*(jj+kk)*g2/lSN/c*lNF)...
                        *besselk(1,sqrt(4*(jj+kk)*g2/lSN/c/lNF));
                    end
                    A4 = A5 - Integral_mu_inf(g2/c/mu_a,1/lNF,(jj+kk)*g2/lSN/c);
                    %
                    A0(ss) = A5;
                    term1 = term1 + (A1 * A2 * A3 * A4);
                end
            end
            term2 = ((1- exp(-mu_a/lSN)).^K)*exp(-g2/lNF/c/mu_a);
            OP_S1_F_ana(aa,rr) =  term1+term2;
        end
    end
end
for xx = 1:length(alpha)
    for yy = 1:length(rho)
        if (0 == isreal(OP_S1_F_ana(xx,yy)))
            OP_S1_F_ana(xx,yy) = 1;
        end
    end
end

%% plot
%
semilogy(rho,OP_S1_F_ana,'*-')
% h = plot(NaN,NaN,'k:',NaN,NaN,'ko');
% legend(h,'Analysis')