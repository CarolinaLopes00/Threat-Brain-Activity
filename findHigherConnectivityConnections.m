%Find the connections with connectivity value higher than a threshold 
%Plot connectivity matrices and brain network scheme

function findHigherConnectivityConnections(condition,methodName,connectivityMatrix,signifChannels,areasOfInterest)

    fprintf('\n\nConnected Channels ('+string(methodName)+', '+string(condition)+'):\n\n')

    connectMatrix=connectivityMatrix;
    if string(methodName)=='Imaginary Part of Coherence' || string(condition)=='Pos-Pre' %to don't account for the connectivity with the itself
        connectMatrix(connectMatrix==0)=[];
    else
        connectMatrix(connectMatrix==1)=[];
    end

    threshold=quantile(connectMatrix,0.65); %quantile 65 for obtaining the threshold for each case 
    
    connectivityValues=[];
    chans=[];

    %find the connectivity values above the value corresponding to the quantile 65%
    for i = 1:length(connectivityMatrix)
        for j = (i+1):length(connectivityMatrix)
            if connectivityMatrix(i,j) > threshold %connectivity value above threshold -> channels are significantly connected
                 connectivityValues=[connectivityValues connectivityMatrix(i,j)];
                 chans=[chans; areasOfInterest(i) areasOfInterest(j)];
            end
        end
    end
       
    auxConnectivityMatrix=connectivityMatrix;
    for i=1:length(auxConnectivityMatrix)
        for j=1:length(auxConnectivityMatrix)
            if i>j
                auxConnectivityMatrix(i,j)=0;
            end
        end
    end

    %order the connectivity values to print the connections in a descendent order
    sortedConnectivityValues=sort(unique(connectivityValues),'descend');
    
    for i=1:length(sortedConnectivityValues)
        [channel1,channel2]=find(auxConnectivityMatrix==sortedConnectivityValues(i));
        fprintf(string(signifChannels(channel1))+'-'+string(signifChannels(channel2))+' ('+string(sortedConnectivityValues(i))+')\n')
    end

    
    %Visualize connectivity matrix (heatmap)
    figure;
    imagesc(connectivityMatrix); 
    colorbar;
    xlabel('Channel'); ylabel('Channel');
    xticklabels(signifChannels); yticklabels(signifChannels);
    title(string(methodName)+' ('+string(condition)+')');

    %Thresholded matrix (for the case Pos-Pre)
    if string(condition)=='Pos-Pre'
        thresholdedMatrix=connectivityMatrix;

        for i = 1:length(connectivityMatrix)
            for j = 1:length(connectivityMatrix)
                if connectivityMatrix(i,j) > threshold %connectivity value above threshold -> channels are significantly connected
                     thresholdedMatrix(i,j)=1;
                else
                    thresholdedMatrix(i,j)=0;
                end
            end
        end

        figure;
        imagesc(thresholdedMatrix); 
        colormap('bone')
        xlabel('Channel'); ylabel('Channel');
        xticklabels(signifChannels); yticklabels(signifChannels);
        title('Thresholded '+string(methodName)+' ('+string(condition)+')');
    end


    % Brain network scheme
    ds.chanPairs=chans;
    ds.connectStrength=connectivityValues;
    ds.connectStrengthLimits=[0 1];
    
    figure; topoplot_connect(ds,'Standard-10-20-Cap19.locs')
    colorbar; title(string(methodName)+' ('+string(condition)+')');

end