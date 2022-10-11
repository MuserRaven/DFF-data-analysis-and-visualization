% DataType = Grating
residual = true; % specify whether you want residual or not 
files = dir('D:\stimulus_data\flashedbars\*.mat');
for f = 10:20  % chooose the number of files you want to extract
   load(files(f).name) 



index = Spike_triggered_average(dataStruct.dendrite.dff,100,30);
% The index of Spine STA when a dendrite Spikes.

%%%%%%%%%%%%%%%%%% CELL47 DENDRITE 2 !!!!
windowsize = 30;
avg = zeros(1,2*windowsize);
storage = nan(length(dataStruct.spine),2*windowsize+1);
%figure

for i = 1:length(dataStruct.spine)
    if residual == true
        dffspine = dataStruct.spine(i).dff_residual;
    else
        dffspine = dataStruct.spine(i).dff_raw;
    end
   
   for j = index(index>=windowsize & length(dffspine)-index(end)>=windowsize)

   
   avg = avg + dffspine(j-windowsize:j+windowsize);
 
   end
   avg = avg./windowsize; storage(i,:) = avg(:,i);
   %plot(avg)
   %hold on 
    
end
%figure
%for i = 1:length(dataStruct.spine)

 %ubplot(length(dataStruct.spine),1,i), imagesc(storage(i,:))
  %xline(windowsize,LineWidth=1.5,Color='r')
 
%end
  
 % xlabel('t')
%  ylabel('dff')
  
figure
imagesc(storage(:,:))
colorbar
xline(30,LineWidth=3)
xlabel('Time')
ylabel('Spine #')
title([files(f).name, ' Spine STA for Dendritic Spike'])

end


