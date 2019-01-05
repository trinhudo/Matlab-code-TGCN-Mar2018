close all
clear all
%
%% Simulation parameters
%
K = 3;                      % # of antenna
rho     = .1:.1:.9;         % power splitting ratio
alpha   = .1:.1:.9;         % time fraction for EH
PS_dB   = 0;                % transmit SNR = Ps/N0 in dB
PS      = 10.^(PS_dB./10);
naN     = (10^(-7))*1e6;    % naN = -100 dBm, BW = 1 MHz
ncN     = (10^(-6))*1e6;    % naN = -90 dBm,  BW = 1 MHz
naF     = (10^(-7))*1e6;
ncF     = (10^(-6))*1e6;
epsilon = 3;                % pathloss exponent
dSF = 10;                   % S-F distance
dSN = 3;
dNF = dSF - dSN;
L   = 1e3;                  % path-loss at reference distance
%
lSN = L*dSN^-3;             % lambda
lSF = L*dSF^-3;
lNF = L*dNF^-3;
%
eta     = 0.7;              % energy conversion coefficient
pN      = 0.1;              % power allocation coefficient
pF      = 1 - pN;
RthN    = 1;                % target data rate of User N bits/s/Hz
RthF    = .1;               % target data rate of User N bits/s/Hz
g2_non  = 2^RthF - 1;
%
SimTimes = 10^5;            % Monte-Carlo repetitions
%
%% Simulation
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
            % random selection
            for yy = 1:SimTimes
                i_rand = randperm(K,1);
                gSsN(yy,1) = gSiN(yy,i_rand);
                gSsF(yy,1) = gSiF(yy,i_rand);
            end
            % Non-cooperative
            snrSsF_non = pF.*PS(ss).*gSsF./(pN.*PS(ss).*gSsF + naF + ncF);
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
            count_coop = 0;
            count_non = 0;
            %
            for zz = 1:SimTimes
                % for non-cooperative
                if (snrSsF_non(zz) < g2_non)
                    count_non = count_non + 1;
                end
                % for Cooperative 
                if (snrSsN_xF(zz) >= g2) && ...
                        (max(snrSsF(zz),snrNF(zz)) < g2)
                    count_coop = count_coop + 1;
                elseif (snrSsN_xF(zz) < g2) && (snrSsF(zz) < g2)
                    count_coop = count_coop + 1;
                end
            end
            OP_coop_random(aa,rr) = count_coop/SimTimes;
            OP_non_random(aa,rr) = count_non/SimTimes;
        end
    end
end

%% plot
surf(alpha,rho,OP_coop_random,'linestyle',':')
hold on
surf(alpha,rho,OP_non_random,'linestyle','-')
%
% zlim([10^-3 10^0])
set(gca, 'ZScale', 'log')
%
xlabel('\rho')
ylabel('\alpha')
zlabel('Outage Probability')
legend('With coop.','Without coop.')
%
set(gca,'XTick',0:.5:1)
set(gca,'YTick',0:.5:1)
