% Kolmogorov-Smirnov test (to find data distribution)

function findDataDistribution(methodName,neutData,stimData,signifChannels)
    
    fprintf("Kolmogorov-Smirnov test ("+string(methodName)+")\n\n") 
    
    for i=1:length(signifChannels)
        for j=(i+1):length(signifChannels)
    
            %neutral data
            h_neut = kstest((squeeze(neutData(i,j,:))-mean(squeeze(neutData(i,j,:))))/std(squeeze(neutData(i,j,:))));  
            if h_neut==1 %h=0 -> normal distribution
                fprintf(signifChannels(i)+"-"+signifChannels(j)+" (neutral) doesn't have a normal distribution\n")
            end
    
            %stimulus data
            h_stim = kstest((squeeze(stimData(i,j,:))-mean(squeeze(stimData(i,j,:))))/std(squeeze(stimData(i,j,:))));
            if h_stim==1 %h=0 -> normal distribution
                fprintf(signifChannels(i)+"-"+signifChannels(j)+" (stimulus) doesn't have a normal distribution\n")
            end
        end
    end

end
