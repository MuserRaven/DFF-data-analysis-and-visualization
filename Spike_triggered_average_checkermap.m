load('checkermap.mat')
windowsize = 10;
net_storage = nan(61,2*windowsize+1);
figure
for f = 1:60 %length(checkermap(1).dff(:,1))  % chooose the number of files you want to extract
   index = Spike_triggered_average(checkermap(2).dff(f,19:end),windowsize,3);
   % The index of Spine STA when a dendrite Spikes.
   avg = zeros(1,2*windowsize+1);
   storage = nan(length(checkermap(2).dff)-18,2*windowsize+1);
   %figure
   
   dffspine = checkermap(2).dff(f,19:17074)';
   
   for j = index(index>=windowsize & length(dffspine)-index(end)>=windowsize)
   avg = avg + dffspine(j-windowsize:j+windowsize);
   end

   avg = avg./windowsize;
   %plot(avg)
   %hold on 
    

%figure
%for i = 1:length(dataStruct.spine)

 %ubplot(length(dataStruct.spine),1,i), imagesc(storage(i,:))
  %xline(windowsize,LineWidth=1.5,Color='r')
 
%end
  
 % xlabel('t')
%  ylabel('dff'
  net_storage(f,:) = avg(:,1)';
end

imagesc(net_storage(:,1:end-1))
colorbar
xline(10,LineWidth=2)
xlabel('Time')
ylabel('Cell')
