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
SimTimes = 5;            % Monte-Carlo repetitions
%
%% Simulation
%
for ss = 1:length(PS_dB)
    disp(strcat('SNR=',num2str(PS_dB(ss)),'dB'));
    for aa = 1:length(alpha)
        disp(strcat('alpha=',num2str(alpha(aa))));
        for rr = 1:length(rho)
            disp(strcat('rho=',num2str(rho(rr))));
            %
            g2 = 2^(RthF*2/(1-alpha(aa))) - 1; % gamma_2
            % channel modelling
            for ii = 1:K
                hiF(:,ii) = sqrt(lSF/2)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
                hiN(:,ii) = sqrt(lSN/2)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            end
            hNF = sqrt(lNF/2)*...
                (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            % channel gains
            giN = abs(hiN.^2);
            giF = abs(hiF.^2);
            gNF  = abs(hNF.^2);
            % find the best far
            [gSsF,I] = max(giF,[],2);
            for yy = 1:SimTimes
                gsN(yy,1) = giN(yy,I(yy));
            end
            % SNR modelling
            snrsNxF = (1-rho(rr)).*pF.*PS(ss).*gsN./...
                ((1-rho(rr)).*pN.*PS(ss).*gsN ...
                + (1-rho(rr))*naN + ncN);
            %
            snrsNxN = (1-rho(rr)).*pN.*PS(ss).*gsN/...
                (1-rho(rr))*naN + ncN;
            %
            snrsF = pF.*PS(ss).*gSsF./(pN.*PS(ss).*gSsF + naF + ncF);
            %
            snrNF = eta.*PS(ss).*gsN.*gNF.*...
                (2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
            % count outage events
            count_2 = 0;
            %
            for zz = 1:SimTimes
                if (snrsNxF(zz) >= g2) && ...
                        (max(snrsF(zz),snrNF(zz)) < g2)
                    count_2 = count_2 + 1;
                elseif (snrsNxF(zz) < g2) && (snrsF(zz) < g2)
                    count_2 = count_2 + 1;
                end
            end
            OP_S2_F_sim(ss) = count_2/SimTimes;
            %% Analytical Results
            a1 = (1-rho(rr))*pF*PS(ss)/((1-rho(rr))*naN + ncN);
            a2 = (1-rho(rr))*pN*PS(ss)/((1-rho(rr))*naN + ncN);
            b1 = pF * PS(ss) / (naF + ncF);
            b2 = pN * PS(ss) / (naF + ncF);
            c  = eta*PS(ss)*(2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
            %
            mu_a = g2/(a1-a2*g2);
            mu_b = g2/(b1-b2*g2);
            %
            c_ = eta*(2*alpha(aa)/(1-alpha(aa))+rho(rr));
            mu_a_ = g2/(1-rho(rr))/(pF-pN*g2);
            SNR = PS(ss)/(naN+ncN);
            % new results
            NewTerm = 0;
            for kk = 0:K
                NewTerm = NewTerm + nchoosek(K,kk)*((-1)^kk)*...
                    exp(-kk*mu_b/lSF) ...
                    *(1 - 1/lSN*Integral_mu_inf(mu_a,1/lSN,g2/lNF/c));
            end
            
            Phi1 = 1 - exp(-g2/lSN/(a1-a2*g2));
            Phi1_asym = mu_a_/lSN/SNR;
            %
            Phi2 = 0;
            %
            for jj = 0:K
                Phi2_temp = nchoosek(K,jj)*((-1)^jj)*...
                    exp(-jj*g2/lSF/(b1-b2*g2));
                Phi2 = Phi2 + Phi2_temp;
            end
            %
            Phi2_asym = (g2/lSF/(pF-pN*g2)/SNR)^K;
            %
            Theta2 = g2/lSN/lNF/c*igamma(0,mu_a/lSN) - ...
                (g2^2)/2/lSN/(lNF^2)/(c^2)/mu_a*exp(-mu_a/lSN) - ...
                (g2^2)/2/(lSN^2)/(lNF^2)/(c^2)*ei(-mu_a/lSN);
            %
            Theta2_asym = -(g2/lSN/lNF/c_/SNR +...
                (g2^2)/2/(lSN^2)/(lNF^2)/(c_^2)/(SNR^2))*...
                (0.5772 + log(mu_a_/lSN/SNR)) -...
                (g2^2)/2/lSN/(lNF^2)/(c_^2)/mu_a_/SNR *...
                (1-mu_a_/lSN/SNR);
            %
            OP_S2_F_ana(ss) = NewTerm;
            OP_S2_F_asym(ss) = Phi2_asym*(Phi1_asym + Theta2_asym);
            
        end
    end
end
for ss = 1:length(PS)
    for aa = 1:length(alpha)
        for rr = 1:length(rho)
            if (0 == isreal(OP_S2_F_ana(ss)))
                OP_S2_F_ana(ss) = 1;
            end
        end
    end
end
%% plot
semilogy(PS_dB,OP_S2_F_sim,'o', ...
    PS_dB,OP_S2_F_ana,'-',...
    PS_dB,OP_S2_F_asym,':')
axis([-20 10 1e-6 1])