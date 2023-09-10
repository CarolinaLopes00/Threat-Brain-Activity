%% -------------------- Preprocessing --------------------

%% Open acquired raw data and remove excess data

%Open and convert EEG file from TRC to mat file
allData = trc_file('C:\Users\lopes\Desktop\aq\1_EEG.TRC'); %load TRC file
allData.get_electrode_info();
[data] = allData.def_data_access(allData.a_n_data_secs, 5, allData.a_file_elec_cell)'; %EEG data

channels = allData.a_file_elec_cell; %channel order
fs = allData.a_samp_freq; %sampling frequency

init_video=228; %instant (in sec) of the begining of the VR video
video_duration=335; %video duration in sec (5:35min= 335 sec)
data=data(init_video*fs:(init_video*fs+(video_duration*fs))-1,:); %remove the instants before and after the video


% Plot raw data
for i=1:length(channels)
    %time
    figure; subplot(2,1,1); 
    plot(0:1/fs:1/fs*(length(data(:,i))-1), data(:,i)); xlabel('Time (s)'); ylabel('Voltage'); title(string(channels(i))+'- Time Domain')
    
    %frequency
    subplot(2,1,2); 
    f = [-numel(data(:,i))/2:numel(data(:,i))/2-1].*fs/numel(data(:,i));
    plot(f, abs(fftshift(fft(data(:,i)))), 'm'); xlabel('Frequency (Hz)'); ylabel('Fourier Transform'); title(string(channels(i))+'- Frequency Domain')
end


% Save the file after removing the instants before and after the video
save('C:\Users\lopes\Desktop\Dados\s1.mat','fs','channels','data') %save in a mat variable


%% Import data (just the data related to the video)

subj=1; %subject number
load('C:\Users\lopes\Desktop\Dados\s'+string(subj)+'.mat','fs','channels','data')

% Plot raw data 
for i=1:length(channels)
    %time
    figure; subplot(2,1,1); 
    plot(0:1/fs:1/fs*(length(data(:,i))-1), data(:,i)); xlabel('Time (s)'); ylabel('Voltage'); title(string(channels(i))+'- Time Domain')
    
    %frequency
    subplot(2,1,2); 
    f = [-numel(data(:,i))/2:numel(data(:,i))/2-1].*fs/numel(data(:,i));
    plot(f, abs(fftshift(fft(data(:,i)))), 'm'); xlabel('Frequency (Hz)'); ylabel('Fourier Transform'); title(string(channels(i))+'- Frequency Domain')
end


%% Filtering

%iirnotch filter
fc=[37 50]; %cut frequencies (remove 37Hz and 50Hz components)
wo=fc/(fs/2); 
bw=wo/35; %bandwidth

[b, a]=iirnotch(wo(1), bw(1)); %notch filter- 37Hz
filtData=filtfilt(b, a, data); %apply filter to the EEG data

[b, a]=iirnotch(wo(2), bw(2)); %notch filter- 50Hz
filtData=filtfilt(b, a, filtData); %apply filter to the EEG data

%bandpass filter 
fc=[0.5 80]; 
filtData=bandpass(filtData, fc, fs, ImpulseResponse="fir"); 


%Plot filtered data 
for i=1:length(channels)
    %time
    figure; subplot(2,1,1); 
    plot(0:1/fs:1/fs*(length(filtData(:,i))-1), filtData(:,i)); xlabel('Time (s)'); ylabel('Voltage'); title(string(channels(i))+'- Time Domain')
    
    %frequency
    subplot(2,1,2); 
    f = [-numel(filtData(:,i))/2:numel(filtData(:,i))/2-1].*fs/numel(filtData(:,i));
    plot(f, abs(fftshift(fft(filtData(:,i)))), 'm'); xlim(fc); xlabel('Frequency (Hz)'); ylabel('Fourier Transform'); title(string(channels(i))+'- Frequency Domain')
end


%% ICA

%Find IC
r=length(channels); %number of independent components to be computed
[Zfica, W, T] = fastICA(filtData',r);


% Plot each IC in the time and frequency domains
for i = 1:r
    figure; subplot(2,1,1); 
    plot(0:1/fs:1/fs*(length(Zfica(i,:))-1), Zfica(i,:),'-'); xlim([0 200]); xlabel('Time (s)'); ylabel('Voltage'); title('Component '+string(i)+'- Time Domain')
    
    subplot(2,1,2); 
    f = [-numel(Zfica(i,:))/2:numel(Zfica(i,:))/2-1].*fs/numel(Zfica(i,:));
    plot(f, abs(fftshift(fft(Zfica(i,:)))), 'm'); xlim(fc); xlabel('Frequency (Hz)'); ylabel('Fourier Transform'); title('Component '+string(i)+'- Frequency Domain')
end

% Plot all ICs in the same graph
figure;
for i = 1:r
    subplot(r,1,i); plot(0:1/fs:1/fs*(length(Zfica(i,:))-1), Zfica(i,:),'-'); xlim([0 50]); ylabel(string(i))
end


% Reconstruction of the EEG data without the noisy IC
c=[3 6 11 16]; % c= noisy components identified in the graphs
T(:,c)=0; %eliminate the noisy components
processedData=(T*Zfica)'; %reconstruct the data without noise


% Plot of reconstructed data after ICA
for i=1:length(channels)
    %time
    figure; subplot(2,1,1); 
    plot(0:1/fs:1/fs*(length(processedData(:,i))-1), processedData(:,i)); xlabel('Time (s)'); ylabel('Voltage'); title(string(channels(i))+'- Time Domain')
    
    %frequency
    subplot(2,1,2); 
    f = [-numel(processedData(:,i))/2:numel(processedData(:,i))/2-1].*fs/numel(processedData(:,i));
    plot(f, abs(fftshift(fft(processedData(:,i)))), 'm'); xlim(fc); xlabel('Frequency (Hz)'); ylabel('Fourier Transform'); title(string(channels(i))+'- Frequency Domain')
end


% Save the preprocessed EEG data
save("C:\Users\lopes\Desktop\Dados\s"+string(subj)+'_processed.mat','fs','channels','processedData') 


%% Plot spectrogram using STFT

for i=1:length(channels) %for all the channels
    figure; spectrogram(processedData(:,i),1/(1/fs),0.5/(1/fs),5/(1/fs),fs,'yaxis'); ylim(fc); title("Short Time Fourier Transform- "+channels(i))
end

