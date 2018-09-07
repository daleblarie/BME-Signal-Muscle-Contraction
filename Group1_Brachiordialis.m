%% EMG tutorial
% Written by Ryan S. McGinnis - ryan.mcginnis@uvm.edu
% August 24, 2018
% EMG signal processing tutorial for BME/EE 296 Biomedical Signal Processing
% Covers basic EMG signal processing and visualization, muscle contration 
% onset detection

%% Load EMG data and Annotations
dt = readtable('elec.csv');
% dt_ann = readtable('Annotations.csv');

%% Convert timestamps to seconds from ms
dt.Timestamp = dt.Timestamp_ms_ / 1000;
g1times_ms = g1times * 1000;
g1times_ms = g1times_ms + dt.Timestamp_ms_(1);
dt.Timestamp = dt.Timestamp - dt.Timestamp(1);


i = 1;
j = 1;
start_ts = zeros(10,1);
while i < length(g1times_ms) 
    start_ts(j) = g1times_ms(i);
    i = i + 2;
    j = j + 1;
end

i = 2;
j = 1;
end_ts = zeros(10,1);
while i < length(g1times_ms) + 1 
    end_ts(j) = g1times_ms(i);
    i = i + 2;
    j = j + 1;
end

i = 1;
start_stop_pairs = [];
while i < length(start_ts)+ 1
    start_stop_pairs = [start_stop_pairs, [start_ts(i); end_ts(i)]];
    i = i + 1;
end
% make a matrix of pairs of start stops, then compare the timestamps to
% start stop points in the matrix and combine the results into one logical
% matrix
% ind_mvc = dt.Timestamp_ms_ >= start_ts & dt.Timestamp_ms_ <= stop_ts;
%% Apply basic EMG signal processing
% Define sampling rate
Fs = 1/mean(diff(dt.Timestamp));

% Apply 4th order high-pass butterworth filter with a cutoff of 20 Hz
n = 2; % filter order (n=2 for a 4th order filter applied with filtfilt)
Wp = 20/(Fs/2); % defining cutoff frequency of 20 Hz in rad/samp
[b,a] = butter(n,Wp,'high'); % From Frigo&Crenna_2009 and Thomas_et_al2014
emg = filtfilt(b,a,dt.Sample_V_); %high-passed signal

% Full wave rectify signal
emg_rect = abs(emg);

% Apply 6th order low-pass elliptic filter with cutoff of 50 Hz
n = 3; % filter order (n=3 for a 6th order filter applied with filtfilt)
Wp = 50/(Fs/2); % defining cutoff frequency of 50 Hz in rad/samp
Rp = 5; % pass band ripple
Rs = 40; % stop band attenuation
[b,a] = ellip(n,Rp,Rs,Wp);
emg_env = filtfilt(b,a,emg_rect);


%% Create Contracting and not contracting part
ind_mvc = zeros(length(dt.Timestamp_ms_),1);
i = 1;
j = 1;
while i < length(start_stop_pairs) + 1
    while j < length(ind_mvc) + 1
        if dt.Timestamp_ms_(j) > start_stop_pairs(2,i) && i < length(start_stop_pairs)
            i = i + 1;
        end
        if ind_mvc(j) ~= 1
            ind_mvc(j) = dt.Timestamp_ms_(j) >= start_stop_pairs(1,i) && dt.Timestamp_ms_(j) <= start_stop_pairs(2,i);
        end    
        j = j + 1;
    end
    i = i + 1;
end
log_mvc = logical(ind_mvc);
%% Visualize enveloped signal
figure; 
plot(dt.Timestamp,emg,'color',[0.6,0.6,0.6]); hold on;
plot(dt.Timestamp,emg_env,'r'); 
legend('raw signal','signal env');
plot(dt.Timestamp,log_mvc*max(emg_env), 'k');

% figure; 
% plot(dt.Timestamp(ind_sig),emg(ind_sig),'color',[0.6,0.6,0.6]); hold on;
% plot(dt.Timestamp(ind_sig),emg_env(ind_sig),'r'); 
% plot(cont_data(:,1),cont_data(:,2)*max(emg(ind_sig)),'k');
% xlabel('Time (s)'); ylabel('Volts'); title('EMG signal and envelop during MVC');
% legend('raw signal','signal env','muscle cont');


%% Determine muscle contractions via Hodges & Bui Method
% Define signal to analyze and time
ind_sig = ind_mvc; % define logical index so MVC data is chosen first
sig = emg_env;
ts = dt.Timestamp;

% Define method parameters
win_len_ms = 25; % window length to consider in ms
num_sd = 1;

% Extract signal-based parameters
emg_env_mvc = emg_env(~log_mvc,:);
baseline = mean(emg_env_mvc(1:round(0.05*Fs),:)); % average of the first 50 ms of data from the mvc trial
sd = std(sig(1:round(0.05*Fs),:)); % standard deviation of the first 50 ms of data from the mvc trial
win_len_samp = round((win_len_ms/1000)*Fs);

% Initialize output array as empty
cont_data = [];

% Determine muscle contraction
stop_ind = win_len_samp;
while stop_ind <= length(sig)
   % Define start of window
   start_ind = stop_ind - win_len_samp + 1;
   
   % Determine if signal window is during contraction
   cont_win = mean(sig(start_ind:stop_ind,:)) > baseline + num_sd*sd; % 1 if muscle contraction, otherwise 0
   
   % Write data to output array
   cont_data = [cont_data; ts(start_ind), cont_win];
   
   % Iterate stop_ind to slide window by 1 sample
   stop_ind = stop_ind + 1;
end
%%  Visualize muscle contraction data
figure; 
plot(dt.Timestamp,emg,'color',[0.6,0.6,0.6]); hold on;
plot(dt.Timestamp,emg_env,'r'); 
plot(cont_data(:,1),cont_data(:,2)*max(emg),'k');
xlabel('Time (s)'); ylabel('Volts'); title('EMG signal and envelop during MVC');
% legend('raw signal','signal env','muscle cont');
