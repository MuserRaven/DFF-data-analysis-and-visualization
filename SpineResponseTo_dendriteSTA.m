% DataType = Grating

index = Spike_triggered_average(dataStruct.dendrite.dff,100,30);
% The index of Spine STA when a dendrite Spikes.

%%%%%%%%%%%%%%%%%% CELL47 DENDRITE 2 !!!!
windowsize = 30;
avg = zeros(1,2*windowsize);
storage = nan(length(dataStruct.spine),2*windowsize+1);
%figure
for i = 1:length(dataStruct.spine)
   dffspine = dataStruct.spine(i).dff_raw;
   for j = index
   if j - windowsize >=0
   avg = avg + dffspine(j-windowsize:j+windowsize);
   else
       disp 'Window Size Too Large'
   end
 
   end
   avg = avg./windowsize; storage(i,:) = avg(:,i);
   %plot(avg)
   %hold on 
    
end
%figure
for i = 1:length(dataStruct.spine)

 subplot(length(dataStruct.spine),1,i), imagesc(storage(i,:))
  xline(windowsize,LineWidth=1.5,Color='r')
 
end
  
  xlabel('t')
  ylabel('dff')
  
figure
imagesc(storage(:,:))
colorbar
xline(30,LineWidth=1.5)
xlabel('Time')
ylabel('Spine #')
title('Cell 47 dendrite2 Spine STA for Dendritic Spike Threshold = 5*std(stim)')


