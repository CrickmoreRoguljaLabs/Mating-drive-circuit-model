% Simulation time (min)
simtime = 15000;

% Time resolution
dt = 1;

% Switches
useNoise = 1;
useBoundary = 1;
useCREB = 1;
useCREBDN = 0;
useStimR = 0;
useStimA = 1;
userecharge = 0;

% Noise level
noisefreq = 1; % /min
noiseamp = 0.15; % in A.U.

% Stimulation properties
stim_ratio = 0.7; % For ratiometric CRN signal (0.6 for 5 min and 0.3 for 1 min)
stim_amp = 4; % For direct subtraction of CRN signal
% stim_times = [4500, 4540, 4580, 4620, 4660, 4700]; % ultraspaced
stim_times = [4500, 4520, 4540, 4560, 4580, 4600]; % Regular
% stim_times = [4500, 4505, 4510, 4515, 4520, 4525]; % Collapsed
% stim_times = 4500; % 1 mating
stim_wind = 5;
stimpoints = ones(stim_wind,1) * stim_times +...
    (0:(stim_wind-1))' * ones(1, length(stim_times));
stimpoints = stimpoints(:);
% stimpoints = 4500:4620;

% NPF properties
NPF_dec = 0.04;
NPF_syn = 0.04;
NPF_act = zeros(simtime,1);

% pCd properties
pCd_dec = 0.02;
pCd_syn = 0.02;
pCd_act = zeros(simtime,1);

% CREB properties
if useCREBDN == 1
    CREB_pow = 0.0001;
else
    CREB_pow = 0.001;
end

CREB_delay = 0;
CREB_conv_pCd = 1;
CREB_conv_NPF = 1;
CREB_act_pCd = zeros(simtime,1);
CREB_act_NPF = zeros(simtime,1);

% Channel properties
Chan_decay = 0.01;
Chan_pow_NPF = 0.002;
Chan_pow_pCd = 0.007;
Chan_delay = 960;
Chan_act_pCd = zeros(simtime,1);
Chan_act_NPF = zeros(simtime,1);

% Artificial recharge properties
rech_amp = 0.05;
rech_wind = [4700,5420]; % Regular
% rech_wind = [4700,4760]; % Short
% rech_wind = [4700,8500]; % Long

% Generate noise
NPF_noisegen = poissonSpikeGen(noisefreq, simtime, 1, 10);
NPF_noise =  mean(reshape(NPF_noisegen,[10,simtime])) * noiseamp;
pCd_noisegen = poissonSpikeGen(noisefreq, simtime, 1, 10);
pCd_noise = mean(reshape(pCd_noisegen,[10,simtime])) * noiseamp;
%
% Track the differentials
dNPF = zeros(simtime, 1);
dpCd = zeros(simtime, 1);
dCREBN = zeros(simtime, 1);
dCREBp = zeros(simtime, 1);
dChanXN = zeros(simtime, 1);
dChanXp = zeros(simtime, 1);

for time = 2 : dt: simtime
    
    % Activity = synaptic excition - activity decay - ChannelX influence
    dNPF(time) = pCd_act(time - 1) * pCd_syn - ...
        NPF_act(time - 1) * NPF_dec - Chan_act_NPF(time - 1) * Chan_pow_NPF;
    dpCd(time) = NPF_act(time - 1) * NPF_syn - ...
        pCd_act(time - 1) * pCd_dec - Chan_act_pCd(time - 1) * Chan_pow_pCd;
    
    % Noise
    if useNoise == 1
        dNPF(time) = dNPF(time) + NPF_noise(time);
        dpCd(time) = dpCd(time) + pCd_noise(time);
    end
        
    % Stim (satiety)
    if sum(time == stimpoints) == 1 && useStimR == 1
        % Ratiometric stimulation
        dNPF(time) = dNPF(time) - NPF_act(time - 1) * (1 - stim_ratio);

    elseif sum(time == stimpoints) == 1 && useStimA == 1
        % Arithmetic stimulation
        dNPF(time) = dNPF(time) - stim_amp;
    end
    
    % Artificial recharge (pcd)
    if time >= rech_wind(1) && time <= rech_wind(2) && userecharge == 1
        dpCd(time) = dpCd(time) + rech_amp;
    end
    
    % CREB update
    if time > CREB_delay && useCREB == 1
        % Adding CREB activity
        dCREBp(time) = dpCd(time) * CREB_conv_pCd;
        dCREBN(time) = dNPF(time) * CREB_conv_NPF;
    end
    
    % Channel X update
    if time > Chan_delay
        dChanXN(time) = CREB_act_NPF(time - Chan_delay) * CREB_pow...
            - Chan_act_NPF(time-1) * Chan_decay;
        dChanXp(time) = CREB_act_pCd(time - Chan_delay) * CREB_pow...
            - Chan_act_pCd(time-1) * Chan_decay;
    end
            
    % Consolidate changes
    NPF_act(time) = NPF_act(time - 1) + dNPF(time) * dt;
    pCd_act(time) = pCd_act(time - 1) + dpCd(time) * dt;
    CREB_act_NPF(time) = CREB_act_NPF(time - 1) + dCREBN(time) * dt;
    CREB_act_pCd(time) = CREB_act_pCd(time - 1) + dCREBp(time) * dt;
    Chan_act_NPF(time) = Chan_act_NPF(time - 1) + dChanXN(time) * dt;
    Chan_act_pCd(time) = Chan_act_pCd(time - 1) + dChanXp(time) * dt;
    
    if useBoundary == 1
        % Limit CREB bounds
        CREB_act_NPF(time) = max(CREB_act_NPF(time), 0);
        CREB_act_pCd(time) = max(CREB_act_pCd(time), 0);
        
        % Limit channel bounds
        Chan_act_NPF(time) = max(Chan_act_NPF(time), 0);
        Chan_act_pCd(time) = max(Chan_act_pCd(time), 0);

        % Actiity floor
        NPF_act(time) = max(NPF_act(time), 0);
        pCd_act(time) = max(pCd_act(time), 0);
    end
end

NPF_act(4800)/NPF_act(4200)

%% Plot
% Calculate time
time = (1 : simtime) / 60;


figure

for i = 1 : 3
    % Subplot
    subplot(3,1,i);
    hold on
    
    % Plot things
    if i == 1
        % Plot activities
        plot(time, NPF_act, time, pCd_act);
        
        % Labels
        legend({'NPF';'pCd'},'Location','SouthEast');
        
    elseif i == 2
        % Plot CREB activities
        if useCREB == 1
            plot(time, CREB_act_NPF, time, CREB_act_pCd)
            legend({'CREB.NPF';'CREB.pCd'},'Location','SouthEast')
        end
        
    elseif i == 3
        % Plot Channel activities
        plot(time, Chan_act_NPF, time, Chan_act_pCd)
        
        % Labels
        legend({'ChanX.NPF';'ChanX.pCd'},'Location','SouthEast')
    end
    
    Ymax = get(gca,'YLIM');
    Ymax = Ymax(2);
    
    % Plot Stim times if needed
    if useStimR == 1
        plot(stim_times/60, Ymax,'.')
    end

    % Plot recharge time if needed
    if userecharge == 1
        plot(rech_wind/60,[Ymax*1.05, Ymax*1.05],'-')
    end

    hold off
    
    ylabel('Activity (A.U.)')
    xlim([0,simtime/60])
    xlabel('Time (hr)')
    xlim([70,160])
end
