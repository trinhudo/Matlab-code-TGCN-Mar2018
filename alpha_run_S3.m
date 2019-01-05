close all
clear all
%
%% Simulation parameters
%
K       = 3;                      % # of antenna
rho     = 0;         % power splitting ratio
alpha   = .1:.01:.9;         % time fraction for EH
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
RthN    = 0.1;                % target data rate of User N bits/s/Hz
RthF    = .1;               % target data rate of User N bits/s/Hz
[pN,pF] = PowerAllocation(RthN,RthF);
%
SimTimes = 10^0;            % Monte-Carlo repetitions
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
            g1 = 2^(2*RthN/(1-alpha(aa))) - 1; % gamma_1
            g2 = 2^(RthF*2/(1-alpha(aa))) - 1; % gamma_2
            % channel modelling
            for ii = 1:K
                hSiF(:,ii) = sqrt(lSF/2)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
                hSiN(:,ii) = sqrt(lSN/2)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            end
            hNF = sqrt(lNF/2)*...
                (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            % channel gains
            gSiN = abs(hSiN.^2);
            gSiF = abs(hSiF.^2);
            gNF  = abs(hNF.^2);
            % Select the best antenna based on the best S-N link
            [gSsN,I] = max(gSiN,[],2);
            for yy = 1:SimTimes
                gSsF(yy,1) = gSiF(yy,I(yy));
            end
            % SNR modelling
            snrSsN_xF = (1-rho(rr)).*pF.*PS(ss).*gSsN./...
                ((1-rho(rr)).*pN.*PS(ss).*gSsN ...
                + (1-rho(rr))*naN + ncN);
            %
            snrSsN_xN = (1-rho(rr)).*pN.*PS(ss).*gSsN/...
                (1-rho(rr))*naN + ncN;
            %
            snrSsF = pF.*PS(ss).*gSsF./(pN.*PS(ss).*gSsF + naF + ncF);
            %
            snrNF = eta.*PS(ss).*gSsN.*gNF.*...
                (2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
            % count outage events
            count_3 = 0;
            %
            for zz = 1:SimTimes
                if (snrSsN_xF(zz) >= g2) && ...
                        (max(snrSsF(zz),snrNF(zz)) < g2)
                    count_3 = count_3 + 1;
                elseif (snrSsN_xF(zz) < g2) && (snrSsF(zz) < g2)
                    count_3 = count_3 + 1;
                end
            end
            OP_S3_F_sim(ss) = count_3/SimTimes;
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
            Phi1 = 0;
            %
            for kk = 0:K
                Phi1_temp = nchoosek(K,kk)*((-1)^kk)*...
                    exp(-kk*mu_a/lSN);
                Phi1 = Phi1 + Phi1_temp;
            end
            %
            Phi1_asym = (mu_a/lSN)^K;
            %
            Phi2 = 1 - exp(-mu_b/lSF);
            %
            Phi2_asym = g2/lSF/(b1-b2*g2);
            %
            Theta2 = 0; 
            Theta2_asym = 0;
            %
            c_tildle = eta*(2*alpha(aa)/(1-alpha(aa))+rho(rr));
            mu_tildle = g2/(1-rho(rr))/(pF-pN*g2);
            SNR = PS(ss)/(naN+ncN);
            %
            for jj = 1:K
                Theta2_temp = nchoosek(K,jj)*((-1)^(jj+1))*...
                    (exp(-jj*mu_a/lSN) - jj/lSN*Integral_mu_inf(mu_a,jj/lSN,g2/lNF/c));
                %
                Theta2_temp_asym = nchoosek(K,jj)*((-1)^(jj+1))*jj*...
                    (g2/lSN/lNF/c*(- 0.5772 - log(jj*mu_a/lSN))-...
                    (g2^2)/2/lSN/(lNF^2)/(c^2)*...
                    ((1-jj*mu_a/lSN)/mu_a +...
                    jj*(0.5772 + log(jj*mu_a/lSN))/lSN));
                %            
                Theta2 = Theta2 + Theta2_temp;
                Theta2_asym = Theta2_asym + Theta2_temp_asym;
            end
            % 
            OP_S3_F_ana(aa,rr) = Phi2*(Phi1 + Theta2); 
        end
    end
end
for ss = 1:length(PS)
    for aa = 1:length(alpha)
        for rr = 1:length(rho)
            if (0 == isreal(OP_S3_F_ana(ss)))
                OP_S3_F_ana(ss) = 1;
            end
        end
    end
end

%% plot
%
semilogy(alpha,OP_S3_F_ana(:,1),'*-')
% h = plot(NaN,NaN,'k:',NaN,NaN,'ko');
% legend(h,'Analysis')