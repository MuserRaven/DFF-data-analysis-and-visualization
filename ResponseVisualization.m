% method 1 linear score = a1*corr + a2*P + error, but at the same time we also need to
% optimize boundary for score. (supervised) 
% method 2 logistic score = 1/exp(-(a1*corr + a2*P + error)) but at the same time we also need to
% optimize boundary for score. (supervised) 
% method 3 k means or knn (unsupervised)
% method 4 random forest or decision tree  (supervised) 

%% Response Visualization And Correlation/P estimate
rng('default')

filelist = dir('D:\stimulus_data\flashedbars\');
%%
figure
for u = 3:20           %length(filelist)
    load(filelist(u).name)


P = nan(length(dataStruct.spine(1).responses),length(dataStruct.spine));
threshold = 1.4; % set the spike detecting threshold
for j = 1:length(dataStruct.spine)

    for i = 1:length(squeeze(dataStruct.spine(1).responses(:,1,1)))
        response = squeeze(dataStruct.spine(j).responses(i,:,:));
        count = max(response,[],2) > threshold;
        count = sum(count)/8; err = std(response,[],1);
        ChatterjeeCorr(response)
        aver = (ChatterjeeCorr(response)+count)/2;
        %subplot(2,8,i); 
        %plot(1:20,response','k',LineWidth=1);ylim([0 6]);axis off; hold on;
        %scatter(1:20,mean(response,1),'b','filled'); hold on; %errorbar(1:20,mean(response,1),err,'r'); % plot the mean of all trials in each timestamp
        %title({['Corr = ', num2str(ChatterjeeCorr(response))],['P-Thres = ',num2str(count)], ...
          %  ['Aver = ',num2str(aver)]})
        scatter(count,ChatterjeeCorr(response),"b")
        hold on 
        
    end
    sgtitle(['Spine Number', num2str(j)])
end
end
% Should I be assuming that the spikes of non-preferred direction are
% uncorrelated ?
% MAYBE I should combine threshold and Chatterjee correlation 
% threshold probability == Number of trials that pass threshold/total number of trials
% If Chatterjee > 0.2 and threshold probability > 0.5 -> consistent success
% If Chatterjee < 0.2 and threshold probability == 0 -> consistent failure
% If Chatterjee < 0.2 and threshold probability > 0 -> inconsistent 
% The above is all I have now for the classification

% But, we can also plot the peaks distribution, std within peaks, k-means, DBSCAN, mean, and all sorts of parameters. 
% 
% But I would like to first try out the threshold + chatterjee to see how well it works!



