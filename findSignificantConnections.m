% Kruskal Wallis test (to find the significant connections)

function findSignificantConnections(methodName,neutData,stimData,signifChannels,avgNeutData,avgStimData,multCorrection)
    
    fprintf('\nSignificant Connections ('+string(methodName)+'):\n\n')
    
    %p-values (from kruskal wallis test)
    p_values = [];
    channPairs=[];
    for i=1:length(signifChannels)
        for j=(i+1):length(signifChannels)
            p_values = [p_values kruskalwallis([squeeze(neutData(i,j,:)) squeeze(stimData(i,j,:))])];
            channPairs=[channPairs; i j];
        end
    end

    %corrected p-values 
    if strcmp(multCorrection,'yes')
        
        %[~, h]=bonf_holm(p_values,0.05); %Holm Bonferroni correction
        [h, ~, ~, ~]=fdr_bh(p_values,0.05,'dep'); %FDR correction
        
        signif_p=find(h);
    
        if length(signif_p) > 0 %if there are significant p values after the correction
            for idx=1:length(signif_p)
                chann_i=channPairs(signif_p(idx),1);
                chann_j=channPairs(signif_p(idx),2);

                if avgStimData(chann_i,chann_j) > avgNeutData(chann_i,chann_j) %we only want the connections with connectivity(stimulus)>connectivity(neutral)
                    fprintf(string(signifChannels(chann_i))+'-'+string(signifChannels(chann_j))+'\n')
                end
            end
        end
  
    else %no multiple comparation correction
        for idx=1:length(p_values)
            chann_i=channPairs(idx,1);
            chann_j=channPairs(idx,2);
            if p_values(idx) < 0.01 && avgStimData(chann_i,chann_j) > avgNeutData(chann_i,chann_j) %significance level=0.01
                fprintf(string(signifChannels(chann_i))+'-'+string(signifChannels(chann_j))+'\n')
            end
        end
    end
end