filelist = dir('C:\Ben Lab\Project1\flashedbars\');

% Sum over all dendrites 
cell = {'cell1_dendrite','cell2_dendrite','cell3_dendrite','cell4_dendrite','cell5_dendrite','cell6_dendrite','cell7_dendrite','cell8_dendrite'};
for c = 3
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
load(filelist(location(1)).name) 
%h = sum(h); %% total number of dendrites in that cell

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
                 threshold = 5*std(dataStruct.spine(i).dff_raw<0);
                 response = squeeze(dataStruct.spine(i).responses(j,:,:)); % each spine's response
                 %net_distance = cat(1,net_distance,dataStruct.spine(i).distance);
                 [M,I] = max(response,[],2);
                 peak(i,1) = mean(M); % value of max
                 peak(i,2) = mean(I); % index of max
                 count = max(response,[],2) > threshold;
                 count = sum(count)/length(dataStruct.spine(1).responses(:,1,1));
                 net_p = cat(1,net_p,count);
                 net_chat = cat(1,net_chat,ChatterjeeCorr(response));           

             end
             net_peak = cat(1,net_peak,peak(:,1));
             net_peak_loc = cat(1,net_peak_loc,peak(:,2));

           end
             subplot(3,length(dataStruct.spine(i).responses(:,1,1)),j)
             scatter(net_peak_loc(net_peak>=threshold),net_peak(net_peak>=threshold),'filled','r')
             hold on 
             scatter(net_peak_loc(net_peak<threshold),net_peak(net_peak<threshold),'filled','b')  
             ylim([0 6])
             ylabel('Response')
             yline(threshold)
             xlabel(['Act = ',num2str(length(net_peak(net_peak>=threshold))/length(net_peak))])
            % hold on 
            % scatter(mean(net_peak_loc(net_peak>=threshold)), length(net_peak(net_peak>=threshold))/length(net_peak),150,'r',"filled")
           %  ylim([0 0.1])
           %  xlim([0 15])
           %  ylabel('% Active')
             xlabel('Mean Peak Loc')

             subplot(3,length(dataStruct.spine(i).responses(:,1,1)),j+length(dataStruct.spine(i).responses(:,1,1)))
             histogram(net_chat,BinWidth=0.05)
             xlabel('Corr(r)')
             xlim([0 1])
             ylim([0 600])
             subplot(3,length(dataStruct.spine(i).responses(:,1,1)),j+2*length(dataStruct.spine(i).responses(:,1,1)))
             histogram(net_p,BinWidth=0.05)
             xlabel('P(r)')
             xlim([0 1])
             ylim([0 600])
            % subplot(4,length(dataStruct.spine(1).responses(:,1,1)),j+3*length(dataStruct.spine(1).responses(:,1,1)))
            % scatter(net_distance,net_p)
            % ylabel('P(r)')
            % xlabel('Distance(um)')

     end
   sgtitle([cellname,' x stimulus'])

end
 
%% SAME AS THE ABOVE BUT 
filelist = dir('C:\Ben Lab\Project1\flashedbars\');

% Sum over all dendrites 

cell = {'cell1_dendrite','cell2_dendrite','cell3_dendrite','cell4_dendrite','cell5_dendrite','cell6_dendrite','cell7_dendrite','cell8_dendrite'
  'cell9_dendrite','cell10_dendrite','cell11_dendrite','cell12_dendrite','cell13_dendrite','cell14_dendrite','cell15_dendrite','cell16_dendrite','cell17_dendrite'
  'cell18_dendrite','cell19_dendrite','cell20_dendrite','cell21_dendrite','cell22_dendrite','cell23_dendrite','cell24_dendrite','cell25_dendrite','cell26_dendrite'
  'cell27_dendrite','cell28_dendrite','cell29_dendrite','cell30_dendrite','cell31_dendrite','cell32_dendrite','cell33_dendrite'};

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
p_spine = [];
corr_spine = [];
peak_loc = [];
peak_spine = [];

for u = location'
  load(filelist(u).name)


for j = 1:length(dataStruct.spine) % number of spine
  threshold = 5*std(dataStruct.spine(j).dff_raw<0);
    peak = nan(length(dataStruct.spine(1).responses(:,1,1)),2);

    for i = 1:length(dataStruct.spine(1).responses(:,1,1)) % Number of stimulus
        disp(length(dataStruct.spine(1).responses(:,1,1)))
        response = squeeze(dataStruct.spine(j).responses(i,:,:));
        count = max(response,[],2) > threshold;
        count = sum(count)/length(dataStruct.spine(1).responses(:,1,1));
        p_spine = cat(1,p_spine,count);
        corr_spine = cat(1,corr_spine,ChatterjeeCorr(response));
        [M,I] = max(response,[],2);
        peak(i,1) = mean(M); % value of max
        peak(i,2) = mean(I); % index of max

    end
    peak_spine = cat(1,peak_spine,peak(:,1));
    peak_loc = cat(1,peak_loc,peak(:,2));
    
end
end

peak_spine = reshape(peak_spine,[length(dataStruct.spine(1).responses(:,1,1)),length(peak_spine)/length(dataStruct.spine(1).responses(:,1,1))]); peak_spine = peak_spine';
peak_loc = reshape(peak_loc,[length(dataStruct.spine(1).responses(:,1,1)),length(peak_loc)/length(dataStruct.spine(1).responses(:,1,1))]); peak_loc = peak_loc';
corr_spine = reshape(corr_spine,[length(dataStruct.spine(1).responses(:,1,1)),length(corr_spine)/length(dataStruct.spine(1).responses(:,1,1))]); corr_spine = corr_spine';
p_spine = reshape(p_spine,[length(dataStruct.spine(1).responses(:,1,1)),length(p_spine)/length(dataStruct.spine(1).responses(:,1,1))]); p_spine = p_spine';
index = nan(length(p_spine(:,1)),1);

for t = 1:length(index)
  if sum(p_spine(t,:)) == 0
    index(t) = 0;
  else
    index(t) = 1;
  end
end
peak_loc = peak_loc(index > 0,:);
peak_spine = peak_spine(index > 0,:);
corr_spine = corr_spine(index > 0,:);
meancorr = mean(corr_spine,1);
p_spine = p_spine(index > 0,:);
meanspine = mean(p_spine,1);
% plot the spine ID vs P(r) Corr(r) average
%plot(1:7,meanspine,'b--o')
%hold on 
%plot(1:7,meancorr,'r--o')
%xlabel('Spine ID')
%legend('Mean P(r)','Mean Corr(r)')
        
%Get the histogram 
for i = 1:length(dataStruct.spine(1).responses(:,1,1))
subplot(3,length(dataStruct.spine(1).responses(:,1,1)),i)
loc = peak_loc(:,i);peak = peak_spine(:,i);
scatter(loc(peak_spine(:,i)>=threshold),peak(peak_spine(:,i)>=threshold),'r','filled');
hold on 
scatter(loc(peak_spine(:,i)<threshold),peak(peak_spine(:,i)<threshold),'b','filled')
xlim([0 20]);ylim([0 6]);yline(threshold);
%xlabel(['%Act = ',num2str(length(loc(peak_spine(:,i)>=threshold))/length(loc))])
subplot(3,length(dataStruct.spine(1).responses(:,1,1)),i+length(dataStruct.spine(1).responses(:,1,1)))
histogram(corr_spine(:,i),BinWidth=0.05)
xlim([0,1])
ylim([0 300])
subplot(3,length(dataStruct.spine(1).responses(:,1,1)),i+2*length(dataStruct.spine(1).responses(:,1,1)))
histogram(p_spine(:,i),BinWidth=0.05)
xlim([0,1])
ylim([0 600])
end
sgtitle([cellname,' x stimulus'])
end
