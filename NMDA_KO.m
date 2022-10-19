load('TSeries-08282022-1146-002.mat')

%% Check Cell's Validity
for j = 1:length(nr_cell)
    threshold = 5*std(nr_cell(j).dff(19:end));
    h = nr_cell(j).dff(19:end);
    if isempty(h(h>=threshold))
        disp 'Invalid Cell'
    else
        disp 'Valid Cell'
    end
end

for j = 1:length(ko_cell)
    threshold = 5*std(ko_cell(j).dff(19:end));
    h = ko_cell(j).dff(19:end);
    if isempty(h(h>=threshold))
        disp 'Invalid Cell'
    else
        disp 'Valid Cell'
    end
end
% cell 31(16 here) NMDA KO is invalid 
%%
figure
subplot(3,1,1)
nr_cell = ce(1:15); ko_cell = ce(16:end);
nr_totaldff = [];
ko_totaldff = [];
for i = 1:length(nr_cell)
    nr_totaldff = cat(1,nr_totaldff,nr_cell(i).dff(19:end));
    ko_totaldff = cat(1,ko_totaldff,ko_cell(i).dff(19:end));
end
histogram(nr_totaldff,100,BinWidth=0.02,EdgeAlpha = 0.5,FaceAlpha = 0.5)
hold on 
histogram(ko_totaldff,100,BinWidth=0.02,EdgeAlpha = 0.5,FaceAlpha = 0.5)
xlim([-0.5 1])
xlabel('Amplitude')
title('Dff Amplitude')
legend('Normal Cell','KO Cell')
%% Frequency Analysis
subplot(3,1,2)
Fs = 30; % sampling frequency    
nr_totaldff = [];

for i = 1:length(nr_cell)
    nr_totaldff = cat(1,nr_totaldff,nr_cell(i).dff(19:end)); 
end
L = length(nr_totaldff);% Length of normal cell signal
ff = fft(nr_totaldff,L);
fffn = ff(1:L/2);
fffn = fffn/max(fffn);
xfft = Fs*(0:L/2-1)/L;
plot(xfft,abs(fffn),Color='red'); hold on

ko_totaldff = [];
for i = 1:length(ko_cell)
    if i ~=31
    ko_totaldff = cat(1,ko_totaldff,ko_cell(i).dff(19:end));
    end
end
L = length(ko_totaldff);% Length of normal cell signal
ff = fft(ko_totaldff,L);
fffk = ff(1:L/2);
fffk = fffk/max(fffk);
xfft = Fs*(0:L/2-1)/L;
plot(xfft,abs(fffk),Color='yellow'); 
xlabel('Normalized Frequency')
ylabel('Normalized Amplitude')
ylim([0,0.1])
legend('Normal','KO')
%%
subplot(3,1,3)
ko_totaldff = [];
for i = 1:length(nr_cell)
    ko_totaldff = cat(1,ko_totaldff,ko_cell(i).dff(19:end));
end

histogram(abs(fffn))
hold on 
histogram(abs(fffk))
legend('Normal Cell','KO Cell')
ylim([0 20000])
xlim([0,0.01])
xlabel('Frequency(Hz)')
title('Dff Frequency')
%% Correlagram
figure
storage = [];
for i = 1:length(ce)
    storage = cat(2,storage,ce(i).dff(19:end));
end
storage = storage';
imagesc(corr(storage'))
yline(15.5,LineWidth=1.5)
xline(15.5,LineWidth=1.5)
colorbar
title('Normal Cell(1-15) vs KO Cells(16-43)')
%% Plot the position 
figure
for i = 1:length(ce)
    if i<=15
      scatter(ce(i).xPos,ce(i).yPos,'filled','b');hold on  % control
    else
      scatter(ce(i).xPos,ce(i).yPos,'filled','r');hold on % treatment
    end
end
