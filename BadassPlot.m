% Example: Cell1 Dendrite 1
net_peak = [];
net_peak_loc = [];

load(filelist(h).name)

for j = 1:length(dataStruct.spine(1).responses(:,1,1)) % The number of single bar stimuli 
   
    subplot(3,length(dataStruct.spine(1).responses(:,1,1)),j)
    peak = nan(length(dataStruct.spine),2);
    for i = 1:length(dataStruct.spine) % the number spine in this dendrite
        response = squeeze(dataStruct.spine(i).responses(j,:,:));
        [M,I] = max(response,[],2);
        peak(i,1) = mean(M); % value of max
        peak(i,2) = mean(I); % index of max
        plot(response','k'); 
        hold on
        scatter(mean(I),mean(M),30,'b','filled');
        hold on 
    end

    net_peak = cat(1,net_peak,mean(peak(:,1)));
    net_peak_loc = cat(1,net_peak_loc,mean(peak(:,2))); 
    scatter(net_peak_loc(j),net_peak(j),100,'g','filled');
    ylim([0,5])
 
 
end
sgtitle('Stimulus x Spine')

chat_total = [];
count_total = [];
for j = 1:length(dataStruct.spine(1).responses(:,1,1)) % The number of single bar stimuli    
  net_count = [];
  net_chatcorr = [];
    for i = 1:length(dataStruct.spine) % the number spine in this dendrite
     
        threshold =  3*std(dataStruct.spine(i).dff_raw<0); % choose conservative threshold
        response = squeeze(dataStruct.spine(i).responses(j,:,:));
        count = max(response,[],2) > threshold;
        count = sum(count)/8;
        net_count = cat(1,net_count,count);

        net_chatcorr = cat(1,net_chatcorr,ChatterjeeCorr(response));
        
        
 
    end
    subplot(3,length(dataStruct.spine(1).responses(:,1,1)),j+length(dataStruct.spine(1).responses(:,1,1)))
    histogram(net_chatcorr,"BinWidth",0.1,'EdgeColor','r','FaceColor',[0.8 0.3 0.8]);
    %scatter(i,net_chatcorr(i),'filled','r'); hold on     
    xlabel('Corr')
    xlim([0 1])
    ylim([0,10])

    subplot(3,length(dataStruct.spine(1).responses(:,1,1)),j+2*length(dataStruct.spine(1).responses(:,1,1)))
    histogram(net_count,8,"BinWidth",0.1,'EdgeColor','r','FaceColor',[0.8 0.3 0.8]);
    %scatter(i,net_count(i),'filled'); hold on 
    xlabel('Probability')
    xlim([0,1])
    ylim([0,10])



   % s = csapi(1:length(dataStruct.spine),net_chatcorr);
   % fnplt(s,'r')


    %subplot(4,length(dataStruct.spine(1).responses(:,1,1)),j+3*length(dataStruct.spine(1).responses(:,1,1)))
    %scatter(net_count,net_chatcorr,10,'filled');
    %xline(0.65,LineWidth=1.5)
    %xline(0.08,LineWidth=1.5)
    %yline(0.144,LineWidth=1.5)
    %xlim([0,1]); ylim([0,1])
    %ylabel('P(r)')
    %xlabel('Corr')

end

%% Sum over all dendrites 

cell = {'cell1_dendrite','cell2_dendrite','cell3_dendrite','cell4_dendrite','cell5_dendrite','cell6_dendrite','cell7_dendrite','cell8_dendrite'};

for c = 8
    figure
cellname = cell(c);
h = [];
for j = 1:527
    if strfind(filelist(j).name,cellname) == 1
       h = cat(1,h,strfind(filelist(j).name,cellname));
    else
       h = cat(1,h,0);
    end
end

location = find(h==1); % where cell1's dendrite locaation 
h = sum(h); %% total number of dendrites in that cell

      for j = 1:length(dataStruct.spine(1).responses(:,1,1)) % The number of single bar stimuli 
              net_chat = [];
              net_p = [];
              net_peak = []; 
              net_peak_loc = [];
              net_distance = [];
           for u = location'
              load(filelist(u).name)          
              peak = nan(length(dataStruct.spine),2);
              
             for i = 1:length(dataStruct.spine) % the number spine in this dendrite
                 treshold = 5*std(dataStruct.spine(i).dff_raw<0);
                 response = squeeze(dataStruct.spine(i).responses(j,:,:)); % each spine's response
                 net_distance = cat(1,net_distance,dataStruct.spine(i).distance);
                 [M,I] = max(response,[],2);
                 peak(i,1) = mean(M); % value of max
                 peak(i,2) = mean(I); % index of max
                 count = max(response,[],2) > threshold;
                 count = sum(count)/length(dataStruct.spine(1).responses(:,1,1)) ;
                 net_p = cat(1,net_p,count);
                 net_chat = cat(1,net_chat,ChatterjeeCorr(response));           

             end
             net_peak = cat(1,net_peak,peak(:,1));
             net_peak_loc = cat(1,net_peak_loc,peak(:,2));

           end
             subplot(3,length(dataStruct.spine(1).responses(:,1,1)),j)
           %  scatter(net_peak_loc(net_peak>=threshold),net_peak(net_peak>=threshold),'filled','r')
           %  hold on 
           %  scatter(net_peak_loc(net_peak<threshold),net_peak(net_peak<threshold),'filled','b')  
          %   ylim([0 6])
          %   ylabel('Response')
            % yline(threshold)
             xlabel(['Act = ',num2str(length(net_peak(net_peak>=threshold))/length(net_peak))])
            % hold on 
             scatter(mean(net_peak_loc(net_peak>=threshold)), length(net_peak(net_peak>=threshold))/length(net_peak),150,'r',"filled")
             ylim([0 0.4])
             xlim([0 15])
             ylabel('% Active')
             xlabel('Mean Peak Loc')

             subplot(3,length(dataStruct.spine(1).responses(:,1,1)),j+length(dataStruct.spine(1).responses(:,1,1)))
             histogram(net_chat,BinWidth=0.05)
             xlabel('Corr(r)')
             xlim([0 1])
             ylim([0 600])
             subplot(3,length(dataStruct.spine(1).responses(:,1,1)),j+2*length(dataStruct.spine(1).responses(:,1,1)))
             histogram(net_p,BinWidth=0.05)
             xlabel('P(r)')
             xlim([0 1])
             ylim([0 600])


     end
   sgtitle([cellname,' x stimulus'])

end
