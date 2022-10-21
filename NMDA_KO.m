load('TSeries-08282022-1146-002.mat')
nr_cell = ce(1:15); ko_cell = ce(16:end);
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
computePeakResp()
% cell 31(16 here) NMDA KO is invalid
%%
figure
x = nr_cell(1).dff(19:end);
threshold = 5*std(x(nr_cell(1).dff(19:end)<0));
t = 1:1:3e2;
f = exp(-t./33);%
plot(deconv(x,f));hold on
yline(threshold)
computePeakResp(nr_cell(1).dff(19:end)')
%
figure
y = deconv(x,f);
plot(y(y>=threshold))




%%
figure
subplot(3,1,1)
nr_totaldff = [];
ko_totaldff = [];
for i = 1:length(nr_cell)
    thres_nr = 5*std(nr_cell(i).dff(19:end));
    thres_ko = 5*std(ko_cell(i).dff(19:end));
    nr_amp = nr_cell(i).dff(19:end);
    ko_amp = ko_cell(i).dff(19:end);
    nr_totaldff = cat(1,nr_totaldff,nr_amp(nr_amp >= thres_nr));
    ko_totaldff = cat(1,ko_totaldff,ko_amp(ko_amp >= thres_ko));
end
histogram(nr_totaldff,100,BinWidth=0.01,EdgeAlpha = 0.5,FaceAlpha = 0.5)
hold on 
histogram(ko_totaldff,100,BinWidth=0.01,EdgeAlpha = 0.5,FaceAlpha = 0.5)
xlim([0 2])
xlabel('Amplitude')
title('Dff Spike Amplitude')
legend(['Normal Cell n = ',num2str(length(nr_totaldff))],['KO Cell n = ',num2str(length(ko_totaldff))])
%% Two-sample Kolmogorov-Smirnov test

[h,p] = kstest2(nr_totaldff,ko_totaldff); 

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
%imagesc(corr(storage'))
histogram(corr(storage'))
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
%% Extract Correlagram clusters and Extract X pos and Y pos and calculate their distance to each other
matrix_corr = corr(storage');
%histogram(matrix_corr)
thres = prctile(matrix_corr,90,'all');
matrix_corr(matrix_corr >= thres) = 1; % get clusters of 1std above the mean
matrix_corr(matrix_corr < thres) = 0;
figure
imagesc(matrix_corr)
title('Correlation Cluster above 90% percentile')

xpos = []; ypos = [];
xpos = cat(1,xpos,ce(1:end).xPos); ypos = cat(1,ypos,ce(1:end).yPos);

for xIndex = 1 : length(xpos)
  for yIndex = 1 : length(ypos)
    matrix_dist(xIndex, yIndex) = sqrt((xpos(xIndex)-xpos(yIndex))^2 + (ypos(xIndex)-ypos(yIndex))^2);
  end
end
thres = prctile(matrix_dist,10,'all');
matrix_dist(matrix_dist < thres ) = 1;
matrix_dist(matrix_dist > thres) = 0;
figure
imagesc(matrix_dist)
title('Distance Cluster below 10% percentile')
%%

overlap = matrix_dist+matrix_corr;
overlap(overlap<2) = 0;
figure

imagesc(overlap)
title('Overlap')


%%
figure
plot(ce(1).dff)
threshold = 5*std(ce(1).dff(19:end));
yline(threshold,LineWidth=1.5)
%%
figure
histogram(ce(1).raw,100)

