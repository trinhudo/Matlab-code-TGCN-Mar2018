close all
clear all
%
%% Simulation parameters
%
K       = 3;                      % # of antenna
rho     = .3;         % power splitting ratio
alpha   = .3;         % time fraction for EH
PS_dB   = 0;                % transmit SNR = Ps/N0 in dB
PS      = 10.^(PS_dB./10);
naN     = (10^(-7))*1e6;    % naN = -100 dBm, BW = 1 MHz
ncN     = (10^(-6))*1e6;    % naN = -90 dBm,  BW = 1 MHz
naF     = (10^(-7))*1e6;
ncF     = (10^(-6))*1e6;
epsilon = 3;                % pathloss exponent
L       = 1e3;                  % path-loss at reference distance
dSF     = 10;                   % S-F distance
dSN     = 1:1:9;
%
eta     = 0.7;              % energy conversion coefficient
RthN    = .1;                % target data rate of User N bits/s/Hz
RthF    = .1;               % target data rate of User N bits/s/Hz
[pN,pF] = PowerAllocation(RthN,RthF);
%
SimTimes = 10^6;            % Monte-Carlo repetitions
%
%% Simulation
%
for ss = 1:length(dSN)
    fprintf('d_{SN} = %d meters \n',dSN(ss))
    dNF     = dSF - dSN(ss);
    %
    lSN     = L*dSN(ss)^-3;             % lambda
    lSF     = L*dSF^-3;
    lNF     = L*dNF^-3;
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
            % SNR modelling
            snrSsN_xF = (1-rho(rr)).*pF.*PS.*gSsN./...
                ((1-rho(rr)).*pN.*PS.*gSsN ...
                + (1-rho(rr))*naN + ncN);
            %
            snrSsN_xN = (1-rho(rr)).*pN.*PS.*gSsN/...
                (1-rho(rr))*naN + ncN;
            %
            snrSsF = pF.*PS.*gSsF./(pN.*PS.*gSsF + naF + ncF);
            %
            snrNF = eta.*PS.*gSsN.*gNF.*...
                (2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
            % count outage events
            count = 0;
            %
            for zz = 1:SimTimes
                %% for DF only
                if (snrSsN_xF(zz) >= g2) && ...
                        (max(snrSsF(zz),snrNF(zz)) < g2)
                    count = count + 1;
                elseif (snrSsN_xF(zz) < g2) && (snrSsF(zz) < g2)
                    count = count + 1;
                end
            end
            OP_S0_random(ss) = count/SimTimes;
        end
    end
end
%% plot
% load('dataSim_far_SNR')
semilogy(dSN,OP_S0_random,'d-')
% hold on
% semilogy(PS_dB,1000./PS.^(K+1))