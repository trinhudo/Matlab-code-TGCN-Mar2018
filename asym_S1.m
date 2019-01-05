close all
clear all
%
%% Simulation parameters
%
K       = 3;                      % # of antenna
rho     = .3;         % power splitting ratio
alpha   = .3;         % time fraction for EH
PS_dB   = -10:5:20;                % transmit SNR = Ps/N0 in dB
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
    fprintf('SNR = %d dB \n',PS_dB(ss))
    for aa = 1:length(alpha)
        for rr = 1:length(rho)
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
            snrFe2e_i(:,:) = min(snriNxF(:,:),max(snrNF(:,:),snriF(:,:)));
            [snrFe2e_b,I] = max(snrFe2e_i,[],2);
            % count outage events
            % method 2
            count = snrFe2e_b < g2;
            OP_S1_F_sim(ss) = sum(count)/SimTimes;
            %% Analysis
            a1 = (1-rho(rr))*pF*PS(ss)/((1-rho(rr))*naN + ncN);
            a2 = (1-rho(rr))*pN*PS(ss)/((1-rho(rr))*naN + ncN);
            b1 = pF * PS(ss) / (naF + ncF);
            b2 = pN * PS(ss) / (naF + ncF);
            c  = eta*PS(ss)*(2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
            mu_a = g2/(a1-a2*g2);
            mu_b = g2/(b1-b2*g2);            
            %
            fun1 = @(y) ((1- exp(-mu_a/lSN)...
                + (1 - exp(-mu_b/lSF)).*...
                (exp(-mu_a/lSN) - exp(-g2/lSN/c./y))).^K)...
                .*1./lNF.*exp(-y./lNF);
            A1 = integral(fun1,0,g2/c/mu_a);
            fun2 = @(y) ((1- exp(-mu_a/lSN)).^K)...
                .*1./lNF.*exp(-y./lNF);
            A2 = integral(fun2,g2/c/mu_a,inf);
            OP_S1_F_ana(ss) =  A1+A2;
            %% Asymptotic
            mu_a_ = g2/(1-rho(rr))/(pF-pN*g2);
            mu_b_ = g2/(pF-pN*g2);
            snr = PS(ss)/(naF+ncF);
            c_ = eta*(2*alpha(aa)/(1-alpha(aa))+rho(rr));
            %
%             % using e^-x = 1 - x
%             fun1 = @(y) (mu_a/lSN + mu_b/lSF - mu_a*mu_b/lSN/lSF ...
%                 - mu_b/lSF.*exp(-g2/lSN/c./y)).^K ...
%                 .* 1/lNF.*exp(-y./lNF); 
%             % separating 1/snr
%             fun1 = @(y) (1/snr*(mu_a_/lSN + mu_b_/lSF) ...
%                 - 1/snr^2*mu_a_*mu_b_/lSN/lSF ...
%                 - 1/snr*mu_b_/lSF.*exp(-g2/lSN/c./y)).^K ...
%                 .* 1/lNF.*exp(-y./lNF);
%             % using the raw integral
%             A1_asym = integral(fun1,0,g2/c/mu_a);
%             % 
            term1 = 0;
            for ii=0:1:K
                for jj=(K-ii):-1:0
                    kk = K - (ii+jj);
                    A1 = factorial(K)/factorial(ii)/factorial(jj)/factorial(kk);
                    A2 = (-1)^(jj+kk);
                    A3 = 1/snr^(ii+2*jj+kk);
                    A4 = (mu_a_/lSN+mu_b_/lSF)^ii;
                    A5 = (mu_a_/lSN)^jj;
                    A6 = (mu_b_/lSF)^(jj+kk);
                    A7 = 1/lNF;
                    %
                    if kk==0
                        A8 = 1;
                    else
%                     % after using 3.324.1
%                     A8 = sqrt(4*kk*g2/lSN/c/(1/lNF))...
%                         *besselk(1,sqrt(4*kk*g2/lSN/c*(1/lNF)));
%                     % separating 1/snr
%                     A8 = sqrt(4*kk*g2/lSN/lNF/c_/snr)*besselk(1,sqrt(4*kk*g2/lSN/lNF/c_/snr));
                    % approximate using xK1(x) = 1+x^2/2ln(2/2)
                    A8 = 1+2*kk*g2/lSN/lNF/c_/snr*log(sqrt(kk*g2/lSN/lNF/c_/snr));
                    end
                    %
%                     % raw form
%                     A9 = A8 - A7*Integral_mu_inf(g2/c/mu_a,1/lNF,kk*g2/lSN/c);
                    % after approximating
                    A9 = A8 - exp(-g2/lNF/c_/mu_a_) ...
                        + kk*g2/lSN/lNF/c_/snr*igamma(0,g2/lNF/c_/mu_a_) ...
                        - kk^2*g2^2/2/lNF/lSN^2/c_^2/snr^2 ...
                        *(c_*mu_a_/g2*exp(-g2/lNF/c_/mu_a_) ...
                        + 1/lNF*ei(-g2/lNF/c_/mu_a_));
                    %
                    term1 = term1 + A1*A2*A3*A4*A5*A6*A9;
                end
            end
            %
            A1_asym = term1;
            A2_asym = 1/(snr^K)*(mu_a_/lSN)^K*exp(-g2/lNF/c_/mu_a_);
            OP_S1_asym(ss) = A1_asym + A2_asym;
            %
        end
    end
end
%% plot
semilogy(PS_dB,OP_S1_F_ana,'k+-',...
    PS_dB,OP_S1_asym,'b--',...
    PS_dB, 1./PS.^(K+1),'r:')
xlabel('SNR (dB)')
ylabel('OP')
legend('Scheme I, (ana.)','Scheme I, (asym.)','(1/SNR)^{K+1}')
% axis([50 100 1e-50 1])