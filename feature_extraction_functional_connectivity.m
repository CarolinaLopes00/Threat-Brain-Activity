%% -------------------- Feature Extraction --------------------

participants=[1 2 6 8 9 10 12 14 16 17 18 19 21 23 24 26 27 29 30 32]; %numbers associated to each participant

%Scared (20 participants): [1 2 6 8 9 10 12 14 16 17 18 19 21 23 24 26 27 29 30 32]
%Not Scared (6 participants): [11 13 15 22 28 31] 

numPartic=length(participants);
numChannels=19;
condition='threat'; %threat, sound or bird


% Feature Extraction

fs=512; %sampling rate

stimulus_length=(round(122.8*fs)):(round(123.8*fs)); %threat stimulus appearance- 122.8s to 123.8s

%stimulusData and neutralData will be posteriorly used in the connectivity analysis
stimulusData=zeros(numChannels,length(stimulus_length),numPartic); %stimulusData variable will contain the EEG time series correspondent to the threat stimulus apperance
neutralData=zeros(numChannels,length(stimulus_length),length((20*fs):fs:(100*fs)),numPartic); %neutralData will contain the EEG time series correspondent to several neutral instants (from 20s to 100s)

bands=["delta","theta","alpha","beta","gamma"]; %bands names
freqBands=[[0.5 4]; [4 7]; [8 12]; [13 30]; [30 80]]; %frequency range for each band

if strcmp(condition,'bird')
    all_subj=zeros(numPartic, 5, 19*length(bands)); %5= number of time instants to be analyzed
else
    all_subj=zeros(numPartic, 4, 19*length(bands)); %4= number of time instants to be analyzed
end


for part=1:numPartic
    load('C:\Users\lopes\Desktop\Dados\s'+string(participants(part))+'_processed.mat','fs','channels','processedData') %load processed data
    
    %saves the temporal data corresponding to the neutral state
    neutIdx=1;
    for neut=(20*fs):fs:(100*fs)
        neutralData(:,:,neutIdx,part)=(processedData(neut:neut+fs,:))'; 
        neutIdx=neutIdx+1;
    end
         
    %saves the temporal data corresponding to the threat stimulus
    stimulusData(:,:,part)=(processedData(stimulus_length,:))'; 

    allChanFeats=[];
    for c=1:length(channels) %for all the channels
        
        [s,f,t,p]=spectrogram(processedData(:,c),1/(1/fs),0.5/(1/fs),5/(1/fs),fs,'yaxis'); %Spectrogram using STFT
        
        %find neutral indices
        [~,neutral_idx1]=min(abs(t-20)); %index corresponding to instant 20 s (start of the neutral part)
        [~,neutral_idx2]=min(abs(t-100)); %index corresponding to instant 100 s (end of the neutral part)

        %find stimulus index
        if strcmp(condition,'threat')
            [~,condit_idx]=min(abs(t-123)); %threat -> 123 s (threat_idx -> [122.5 123.5]s)
        elseif strcmp(condition,'sound')
             [~,condit_idx]=min(abs(t-194)); %sound -> 194 s (sound_idx -> [193.5 194.5]s)
        else %bird
            [~,condit_idx]=min(abs(t-286)); %bird -> 285 s (bird_idx -> [284.5 285.5]s)
        end
       
        %Indices corresponding to the time segments of interest
        segments=[neutral_idx1 condit_idx condit_idx+1 condit_idx+2]; %1 before stimulus and 3 after stimulus
        if strcmp(condition,'bird')
            segments=[neutral_idx1 condit_idx condit_idx+1 condit_idx+2 condit_idx+3]; %1 before stimulus and 4 after stimulus
        end

        featsArray=[];
        for i=1:length(segments) %for each time segment
            
            feats=zeros(1,length(freqBands));
            
            for b=1:length(freqBands) %for each frequency band
                if i==1 %i=1 corresponds to the neutral (before stimulus) time segment
                    for ii=neutral_idx1:neutral_idx2

                        %absolute power of the considered frequency band 
                        feats(b)=feats(b)+sum(p([find(f >= freqBands(b,1) & f <= freqBands(b,2))],ii));             
                    end

                    %mean of (neutral_idx1:neutral_idx2) neutral segments
                    feats(b)=feats(b)/(length(neutral_idx1:neutral_idx2)); 
                end

                %absolute power of the considered frequency band    
                feats(b)=sum(p([find(f >= freqBands(b,1) & f <= freqBands(b,2))],segments(i)));          
            end

            %relative power of each frequency band
            feats=feats/sum(p([find(f >= 0.5 & f <= 80)],segments(i))); 
            featsArray=[featsArray; feats];
        end
        allChanFeats=[allChanFeats featsArray];
    end

    all_subj(part,:,:)=allChanFeats;

    fprintf(string(part))
end

columnNames=[]; %column names for the final features table
for c=1:length(channels)
    for b=1:length(bands)
        columnNames=[columnNames channels(c)+"_"+bands(b)];
    end
end

% create table with column names
all_subj_table=table(all_subj);
all_subj_feats=all_subj_table.Variables;


%% Boxplots for each frequency band and channel, englobing all the participants

for i=1:length(channels)
    for j=1:length(bands) 
        label=channels(i)+"_"+bands(j);
        figure; boxplot([all_subj_feats(:,1,find(columnNames==label)), all_subj_feats(:,2,find(columnNames==label)), all_subj_feats(:,3,find(columnNames==label)), all_subj_feats(:,4,find(columnNames==label))],'symbol', '','Labels',{'t(neutral)','t(stimulus)','t(stimulus+1)','t(stimulus+2)'})
        title("Channel: "+channels(i)+"; Frequency Band: "+bands(j)); ylabel('Relative PSD'); 
        
    end
end


%% Kolmogorov-Smirnov test (to find data distribution)

fprintf("Kolmogorov-Smirnov test\n\n") 

for i=1:length(channels)
    for j=1:length(bands) 
        for k=1:length(segments)
            label=channels(i)+"_"+bands(j);
            
            data=all_subj_feats(:,k,find(columnNames==label));
            h = kstest((data-mean(data))/std(data));  
            
            if h==1 %h=0 -> normal distribution
                fprintf("Group "+channels(i)+"_"+bands(j)+" (segment"+string(k)+") doesn't have a normal distribution\n")
            end
        end
    end
end


%% Kruscal-Wallis test (find relevant brain areas and frequency bands)

fprintf("Relevant features based on Kruscal-Wallis\n\n")

for i=1:length(channels)
    for j=1:length(bands) 
        label=channels(i)+"_"+bands(j);

        figure; 

        if strcmp(condition,'bird')
            [p,tbl,stats] = kruskalwallis([all_subj_feats(:,1,find(columnNames==label)) all_subj_feats(:,2,find(columnNames==label)) all_subj_feats(:,3,find(columnNames==label)) all_subj_feats(:,4,find(columnNames==label)) all_subj_feats(:,5,find(columnNames==label))], {'t(neutral)','t(stimulus)','t(stimulus+1)','t(stimulus+2)','t(stimulus+3)'},'off');
        else
            [p,tbl,stats] = kruskalwallis([all_subj_feats(:,1,find(columnNames==label)) all_subj_feats(:,2,find(columnNames==label)) all_subj_feats(:,3,find(columnNames==label)) all_subj_feats(:,4,find(columnNames==label))], {'t(neutral)','t(stimulus)','t(stimulus+1)','t(stimulus+2)'},'off');
        end
        
        %multiple comparation (with Bonferroni correction)
        [c,~,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni"); title("Channel: "+channels(i)+"; Frequency Band: "+bands(j)); 

        idx=find(c(:,6)<=0.05); 
        %c(:,6)-> column with the p values
        %c(:,6)<=0.05-> p values <= significance level 0.05

        if length(idx)>0 %if there are p values < 0.05
            for k=1:length(idx)
                if c(idx(k),1)==1 || c(idx(k),2)==1 %because we want the significant difference compared with the instant before stimulus (neutral condition)
                    fprintf(label+" ("+string(c(idx(k),1))+","+string(c(idx(k),2))+")\n")
                end
            end
        end
    end
end


%% Significative channels 

%Areas of interest (channels) obtained from the Kruscal-Wallis test
signifChannels={'Fp1','Fp2','F7','Fz','F8','T3','C3','Cz','T4','T5'}; %scared

%Indices of brain areas of interest (channels)
areasOfInterest=[];
for i=1:length(signifChannels)
    areasOfInterest = [areasOfInterest find(string(channels)==string(signifChannels(i)))];
end

numRelevChann=length(areasOfInterest); %number of areas/channels of interest (discriminative between the neutral and stimulus conditions)

%Remove the irrelevant channels 
neutralData = neutralData(areasOfInterest,:,:,:);
stimulusData = stimulusData(areasOfInterest,:,:);




%% -------------------- Functional Connectivity --------------------

%Frequency band to be analysed 
bandOfInterest='all'; %all, delta, theta, alpha, beta, gamma
freqRange=[0.5 80]; %[0.5 80], [0.5 4], [4 7], [8 12], [13 30], [30 80]


%% Coherence and Imaginary Part of Coherence

%Pre Stimulus (neutral condition)
preCoherence = zeros(numRelevChann, numRelevChann,size(neutralData,3), numPartic);
preImagCoherence = zeros(numRelevChann, numRelevChann,size(neutralData,3), numPartic);

%Pos Stimulus (threat condition)
posCoherence = zeros(numRelevChann, numRelevChann, numPartic);
posImagCoherence = zeros(numRelevChann, numRelevChann, numPartic);

for p = 1:numPartic
    for i = 1:numRelevChann
        for j = 1:numRelevChann
            
            %Pre Stimulus
            for neutIdx=1:size(neutralData,3)

                %cross-spectral power density
                [pre_Sxy f] = cpsd(squeeze(neutralData(i,:,neutIdx,p)), squeeze(neutralData(j,:,neutIdx,p)),[],[],2048,fs); 
                [pre_Sxx f] = cpsd(squeeze(neutralData(i,:,neutIdx,p)), squeeze(neutralData(i,:,neutIdx,p)),[],[],2048,fs); 
                [pre_Syy f] = cpsd(squeeze(neutralData(j,:,neutIdx,p)), squeeze(neutralData(j,:,neutIdx,p)),[],[],2048,fs); 
                       
                %restrict to the frequency band of interest
                freqIdx = find(f >= freqRange(1) & f <= freqRange(2));
                pre_Sxy=pre_Sxy(freqIdx);
                pre_Sxx=pre_Sxx(freqIdx);
                pre_Syy=pre_Syy(freqIdx);
                
                preCoherency = pre_Sxy ./ sqrt(pre_Sxx .* pre_Syy); %Coherency (normalized cross-spectrum) 
                preCoherence(i, j, neutIdx, p) = mean(abs(preCoherency)); %Coherence (absolute value of coherency)
                preImagCoherence(i, j, neutIdx, p) = mean(abs(imag(preCoherency))); %Imaginary part of coherence
            end
            
            %Pos Stimulus

            %cross-spectral power density
            [pos_Sxy f] = cpsd(squeeze(stimulusData(i,:,p)), squeeze(stimulusData(j,:,p)),[],[],2048,fs); 
            [pos_Sxx f] = cpsd(squeeze(stimulusData(i,:,p)), squeeze(stimulusData(i,:,p)),[],[],2048,fs); 
            [pos_Syy f] = cpsd(squeeze(stimulusData(j,:,p)), squeeze(stimulusData(j,:,p)),[],[],2048,fs); 
                  
            %restrict to the frequency band of interest
            pos_Sxy=pos_Sxy(freqIdx);
            pos_Sxx=pos_Sxx(freqIdx);
            pos_Syy=pos_Syy(freqIdx);
            
            posCoherency = pos_Sxy ./ sqrt(pos_Sxx .* pos_Syy);
            posCoherence(i, j, p) = mean(abs(posCoherency));
            posImagCoherence(i, j, p) = mean(abs(imag(posCoherency))); 
        end
    end
end

% Average of connectivity across the several neutral instants
preCoherence=squeeze(mean(preCoherence,3));
preImagCoherence=squeeze(mean(preImagCoherence,3));

% Average of all the participants
avgPreCoherence=mean(preCoherence,3);
avgPreImagCoherence=mean(preImagCoherence,3);
avgPosCoherence=mean(posCoherence,3);
avgPosImagCoherence=mean(posImagCoherence,3);

% Pos-Pre (difference between the connectivity values of the neutral and threat conditions)
diffCoherence=avgPosCoherence-avgPreCoherence; 
diffImagCoherence=avgPosImagCoherence-avgPreImagCoherence; 


%% kstest to find data distribution
findDataDistribution('Coherence',preCoherence,posCoherence,signifChannels)
fprintf("-----\n")
findDataDistribution('Imaginary Part of Coherence',preImagCoherence,posImagCoherence,signifChannels)


%% find connections that present significant differences between the neutral and threat conditions
findSignificantConnections('Coherence',preCoherence,posCoherence,signifChannels,avgPreCoherence,avgPosCoherence,'yes')
findSignificantConnections('Imaginary Part of Coherence',preImagCoherence,posImagCoherence,signifChannels,avgPreImagCoherence,avgPosImagCoherence,'yes')


%% find channels pairs with higher difference values between the connectivity of the neutral and threat conditions

%Pre (COH)
findHigherConnectivityConnections('Pre Stimulus','Coherence',avgPreCoherence,signifChannels,areasOfInterest)
%Pre (ICOH)
findHigherConnectivityConnections('Pre Stimulus','Imaginary Part of Coherence',avgPreImagCoherence,signifChannels,areasOfInterest)

%Pos (COH)
findHigherConnectivityConnections('Pos Stimulus','Coherence',avgPosCoherence,signifChannels,areasOfInterest)
%Pos (ICOH)
findHigherConnectivityConnections('Pos Stimulus','Imaginary Part of Coherence',avgPosImagCoherence,signifChannels,areasOfInterest)

%Pos-Pre (COH)
findHigherConnectivityConnections('Pos-Pre','Coherence',diffCoherence,signifChannels,areasOfInterest)
%Pos-Pre (ICOH)
findHigherConnectivityConnections('Pos-Pre','Imaginary Part of Coherence',diffImagCoherence,signifChannels,areasOfInterest)


%% Weighted Phase Lag Index

%Pre stimulus (neutral condition)
pre_wPLI = ones(numRelevChann, numRelevChann, size(neutralData,3), numPartic);

%Pos stimulus (threat condition)
pos_wPLI = ones(numRelevChann, numRelevChann, numPartic);


for p=1:numPartic
    for i = 1:numRelevChann
        for j = 1:numRelevChann
            if i~=j

                %Pre Stimulus 
                for neutIdx=1:size(neutralData,3)

                    %cross-spectral power density
                    [pre_cpsd f] = cpsd(squeeze(neutralData(i,:,neutIdx,p)), squeeze(neutralData(j,:,neutIdx,p)),[],[],2048,fs); 
                    
                    %restrict to the frequency band of interest
                    freqIdx = find(f >= freqRange(1) & f <= freqRange(2));
                    pre_cpsd=pre_cpsd(freqIdx);
                    
                    pre_wPLI(i, j, neutIdx, p) = abs((sum(abs(imag(pre_cpsd)).*sign(imag(pre_cpsd)))))/(sum(abs(imag(pre_cpsd))));   
                end

                %Pos Stimulus
                pos_cpsd = cpsd(squeeze(stimulusData(i,:,p)), squeeze(stimulusData(j,:,p)),[],[],2048,fs);
                
                freqIdx = find(f >= freqRange(1) & f <= freqRange(2));
                pos_cpsd=pos_cpsd(freqIdx);
                
                pos_wPLI(i, j, p) = abs((sum(abs(imag(pos_cpsd)).*sign(imag(pos_cpsd)))))/(sum(abs(imag(pos_cpsd))));   
            end
        end
    end
end

%average of all neutral instants
pre_wPLI=squeeze(mean(pre_wPLI,3)); 

%average of all the participants
avgPre_wPLI=mean(pre_wPLI,3); 
avgPos_wPLI=mean(pos_wPLI,3);

%Pos-Pre (difference between the connectivity values of the neutral and threat conditions)
diffwPLI=avgPos_wPLI-avgPre_wPLI; 


%% kstest to find data distribution
findDataDistribution('Weighted Phase Lag Index',pre_wPLI,pos_wPLI,signifChannels)


%% find connections that present significant differences between the neutral and threat conditions
findSignificantConnections('Weighted Phase Lag Index',pre_wPLI,pos_wPLI,signifChannels,avgPre_wPLI,avgPos_wPLI,'yes')


%% find channels pairs with higher difference values between the connectivity of the neutral and threat conditions

%Pre (wPLI)
findHigherConnectivityConnections('Pre Stimulus','Weighted Phase Lag Index',avgPre_wPLI,signifChannels,areasOfInterest)

%Pos (wPLI)
findHigherConnectivityConnections('Pos Stimulus','Weighted Phase Lag Index',avgPos_wPLI,signifChannels,areasOfInterest)

%Pos-Pre (wPLI)
findHigherConnectivityConnections('Pos-Pre','Weighted Phase Lag Index',diffwPLI,signifChannels,areasOfInterest)


%% Mean Phase Coherence

%Pre stimulus (neutral condition)
pre_mpc = zeros(numRelevChann, numRelevChann, size(neutralData,3), numPartic);

%Pos stimulus (threat condition)
pos_mpc = zeros(numRelevChann, numRelevChann, numPartic);


for p = 1:numPartic
    for i = 1:numRelevChann
        for j = 1:numRelevChann
            
            %Pre Stimulus
            for neutIdx=1:size(neutralData,3)

                %Phase for the two channels
                pre_phase1=angle(hilbert(squeeze(neutralData(i,:,neutIdx,p)))); 
                pre_phase2=angle(hilbert(squeeze(neutralData(j,:,neutIdx,p))));

                %Phase differences between the two channels
                pre_phaseDiff = pre_phase1 - pre_phase2; 
        
                pre_mpc(i,j,neutIdx,p) = abs(mean(exp(1i*pre_phaseDiff)));
            end

            %Pos Stimulus

            %Phase for the two channels
            pos_phase1=angle(hilbert(squeeze(stimulusData(i, :, p)))); 
            pos_phase2=angle(hilbert(squeeze(stimulusData(j, :, p))));

            %Phase differences between the two channels
            pos_phaseDiff = pos_phase1 - pos_phase2;
    
            pos_mpc(i,j,p) = abs(mean(exp(1i*pos_phaseDiff)));
        end
    end
end

%Average of all the neutral instants
pre_mpc=squeeze(mean(pre_mpc,3)); 

%Average the MCP values across participants
avgPre_mpc = mean(pre_mpc, 3);
avgPos_mpc = mean(pos_mpc, 3);

%Pos-Pre (difference between the connectivity values of the neutral and threat conditions)
diffMPC=avgPos_mpc-avgPre_mpc; 


%% kstest to find data distribution
findDataDistribution('Mean Phase Coherence',pre_mpc,pos_mpc,signifChannels)


%% find connections that present significant differences between the neutral and threat conditions
findSignificantConnections('Mean Phase Coherence',pre_mpc,pos_mpc,signifChannels,avgPre_mpc,avgPos_mpc,'yes')


%% find channels pairs with higher difference values between the connectivity of the neutral and threat conditions

%Pre (MPC)
findHigherConnectivityConnections('Pre Stimulus','Mean Phase Coherence',avgPre_mpc,signifChannels,areasOfInterest)

%Pos (MPC)
findHigherConnectivityConnections('Pos Stimulus','Mean Phase Coherence',avgPos_mpc,signifChannels,areasOfInterest)

%Pos-Pre (MPC)
findHigherConnectivityConnections('Pos-Pre','Mean Phase Coherence',diffMPC,signifChannels,areasOfInterest)


%% Directed Transfer Function (eConnectome)
%NOTE:this connectivity method was implemented using the eConnectome toolbox

for part=1:numPartic
    load('C:\Users\lopes\Desktop\Dados\s'+string(participants(part))+'_processed.mat','processedData') %load processed data of each participant
    stimData=(processedData(stimulus_length,areasOfInterest))'; %temporal data corresponding to the threat stimulus

    EEG.data=stimData;
    EEG.labels=signifChannels;
    EEG.type='EEG';
    EEG.nbchan=numRelevChann;
    EEG.points=length(stimData);
    EEG.srate=fs;
    EEG.labeltype='standard';
    EEG.unit='V';
    
    save('C:\Users\lopes\Desktop\DataStimulus\stim'+string(part)+'.mat','EEG') %save in a mat variable compatible with the eConnectome toolbox
end


%% Average DTF matrix across participants

avgDTF=zeros(numRelevChann,numRelevChann,numPartic);

for part=1:numPartic
    load('C:\Users\lopes\Desktop\DataDTF\dtf'+string(part)+'.mat','DTF') %load dtf file obtained from eConnectome
    avgDTF(:,:,part)=mean(DTF.matrix,3);
end

avgDTF=mean(avgDTF,3);

%Visualize connectivity matrix (heatmap)
figure;
imagesc(avgDTF); 
colorbar;
xlabel('Channel'); ylabel('Channel');
xticklabels(signifChannels); yticklabels(signifChannels);
title('Directed Transfer Function (Pos Stimulus)');

