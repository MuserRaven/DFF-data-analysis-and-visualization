
% Example: Cell1 Dendrite 2
figure
net_peak = [];
net_peak_loc = [];
for j = 1:length(dataStruct.spine(1).responses(:,1,1)) % The number of single bar stimuli 
    
    subplot(3,length(dataStruct.spine(1).responses(:,1,1)),j)
    peak = nan(length(dataStruct.spine),2);
    for i = 1:length(dataStruct.spine) % the number spine in this dendrite

        response = squeeze(dataStruct.spine(i).responses(j,:,:));
        [M,I] = max(response,[],2);
        peak(i,1) = mean(M); % value of max
        peak(i,2) = mean(I); % index of max
        plot(response','k'); %axis off;
        hold on
        scatter(mean(I),mean(M),200,'b','filled')
        hold on 
    end
    net_peak = cat(1,net_peak,mean(peak(:,1)));
    net_peak_loc = cat(1,net_peak_loc,mean(peak(:,2)));
    scatter(net_peak_loc(j),net_peak(j),500,'g','filled')
end
sgtitle('Stimulus x Spine')

%%

for j = 1:length(dataStruct.spine(1).responses(:,1,1)) % The number of single bar stimuli    
 net_count = [];
    for i = 1:length(dataStruct.spine) % the number spine in this dendrite
        response = squeeze(dataStruct.spine(i).responses(j,:,:));
        count = max(response,[],2) > threshold;
        count = sum(count)/8;
        net_count = cat(1,net_count,count);
        subplot(3,length(dataStruct.spine(1).responses(:,1,1)),j+length(dataStruct.spine(1).responses(:,1,1)))
        %histogram(net_count,8,"BinWidth",0.1,'EdgeColor','r','FaceColor',[0.8 0.3 0.8]);
        plot(net_count,'o'); hold on 
        %ylabel('Spine Count')
        %xlabel('Probability')
        %xlim([0,1])
        %ylim([0,22])
        ylabel('P(r)')
        xlabel('Spine ID')
        xlim([0,length(dataStruct.spine)])
        ylim([0,1])
        
    end
    hold on 
    pp = csapi(1:length(dataStruct.spine),net_count);
    fnplt(pp)
    legend('data','cubic spline')
end

%%
for j = 1:length(dataStruct.spine(1).responses(:,1,1)) % The number of single bar stimuli    
 net_chatcorr = [];
    for i = 1:length(dataStruct.spine) % the number spine in this dendrite
        response = squeeze(dataStruct.spine(i).responses(j,:,:));
        count = max(response,[],2) > threshold;
        count = sum(count)/8;
        net_chatcorr = cat(1,net_chatcorr,ChatterjeeCorr(response));
        subplot(3,length(dataStruct.spine(1).responses(:,1,1)),j+2*length(dataStruct.spine(1).responses(:,1,1)))
        %histogram(net_chatcorr,"BinWidth",0.1,'EdgeColor','r','FaceColor',[0.8 0.3 0.8]);
        plot(net_chatcorr,'o'); hold on 
        %ylabel('Spine Count')
        %xlabel('Corr')
        %xlim([0,1])
        %ylim([0,22])
        ylabel('Corr(r)')
        xlabel('Spine ID')
        xlim([0,length(dataStruct.spine)])
        ylim([0,1])
        
    end
    hold on 
    s = csapi(1:length(dataStruct.spine),net_chatcorr);
    fnplt(s)
    legend('data','cubic spline')
end
