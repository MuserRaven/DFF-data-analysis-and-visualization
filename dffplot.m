
%%%% This page is for the data analysis for dff data



%% STORE DFF INFO IN A MATRIX DFFMATRIX
Numberofcell = length(ce);
dffmatrix = nan(length(ce(1).dff),Numberofcell); 
% matrix that store dff infor of all neurons in ROI over all timesteps
meandff = nan(length(ce),1);
% mean of dff for each neurons in ROI
for i = 1:Numberofcell
     dffmatrix(:,i)=ce(i).dff;    
end
dffmatrix = rmmissing(dffmatrix);
Timelaspe = length(dff)-length(dffmatrix);
% the timelaspse for the neurons to respond

for i = 1:Numberofcell
    meandff(i)=mean(dffmatrix(:,i)); % get the mean dff for all neurons
end
plot(meandff)
xlabel('# of Neurons')
ylabel('Mean DFF')
figure
cell = 9; % the number of neuron we are interested in
plot(dffmatrix(:,cell))
movavgvector = ones(100,1); % number of data point used for each bin 
num = (1/length(movavgvector))*movavgvector; % the moving average filter
den = 1;
y = filter(num,den,dffmatrix(:,cell)); % apply a moving average filter to create a baseline 
hold on 
plot(y,'r')
xlabel('Time Steps')
ylabel('DFF and Baseline')
legend('DFF','MA Baseline')

% do some fourier transform
figure,plot(autocorr(dffmatrix(:,cell),NumLags=100))
% do some auto correlation
xlabel('Time Lag')
ylabel('AutoCorrelation')

filt = filter2(den,dffmatrix(:,1));
plot(filt)



%% Processed Df/f data
title('Single-Sided Amplitude Spectrum of X(t)')

subplot(211),plot(dff)
xlabel('Time Series')
ylabel('DF(t)-F(0)/F(0)')
normdff = normalize(dff);
subplot(212), plot(normdff)
xlabel('Time Series')
ylabel('DF(t)-F(0)/F(0)')

dff = rmmissing(dff);
fdff = fft(dff); % fourier transform
scatterplot(fdff)
xlabel('Frequency(Hz)')
ylabel('Amplitude')
xlim([-100,100])
ylim([-100,100])

figure,plot(autocorr(dff,NumLags=100))
xlabel('Time Lag')
ylabel('Correlation R')

disp('number of sample scan:')
disp(size(dj));

%meandff = mean(dff); doing a gaussian curve fit to the histgram of dff
 a = 0.3;%  standard deviation
 b = 0; %  mean
 ny = length(dff);
 y = dff;   % data
 mu = mean(y);    % data mean
 sd = std(y);     % data std
 nbins = round(ny/100);
 hg = histogram(y, nbins);
 hold on;
 y_bin = hg.Values;
 x_bin = (hg.BinEdges(1:end-1)+hg.BinEdges(2:end))/2;
 y_pdf = 1/(2*pi*sd)*exp(-(x_bin-mu).^2/(2*sd^2));
 area_hist = trapz(x_bin, y_bin);
 area_pdf = trapz(x_bin, y_pdf);
 h_fit = plot(x_bin,y_pdf*area_hist/area_pdf,'LineWidth',3);
 legend(h_fit, sprintf('mu %.3f, std %.3f', mu, sd));
 title(sprintf('Gaussian fit to histogram of %d observations with %d bins', length(y), nbins));
%% Raw Voltage Data
fid = readtable('TSeries-08282022-1146-002_Cycle00001_VoltageRecording_001.csv'); %read the voltage file
voltage = fid{:,:}; % convert table to matrix
ioi = 1:200; % interval of interest
num = (1/length(ones(5,1)))*ones(5,1);

%subplot(221),plot(voltage(ioi,1),voltage(ioi,3)),hold on, plot(filter(num,den,voltage(1:length(ioi)/(2*length(num)),3)))
%subplot(222),plot(fft(voltage(1:1000,3)))

%% Correlation Matrix for Different Cells' Normalized DFF
target = dffmatrix;
corr = corrcoef(target);
figure
imagesc(corr) 
colorbar
title('Correlation Matrix of All Cells in ROI')
xlabel('cell #')
ylabel('cell #')

%% Auto/Cross Correlation for cells
ncell = 5; % Number of Cell
NumberLag = 50; % Number of Lag
figure 
for i = 1:ncell
   subplot(ncell,1,i)
   plot(autocorr(dffmatrix(:,i),NumLags=NumberLag))
end
figure
for i = 1:ncell
    subplot(ncell,1,i)
    [c,lag] = xcorr(dffmatrix(:,i),(dffmatrix(:,i+1)));
    plot(lag,c)
end
%% Visualization of The Entire dffmatrix 
figure
imagesc(dffmatrix(1:100,:))
colorbar
%% Median and P-percentile mask and Moving Average on Raw Data(Reconstructing the getTraces Dff ) 
K_neighbor = 5; % neibor size for mask
percent = 0.08; % percentile for mask
n = 1; % the label of the neuron 
empty_mask = zeros(K_neighbor,1);
raw_vector = [ce(n).raw];
raw_vector = raw_vector(31:end);
percent_mask = nan(length(raw_vector)-2*K_neighbor,1); % to remove unusable data on edges
for p = K_neighbor+1:(length(raw_vector)-(K_neighbor+1))
percent_mask(p) =  prctile(raw_vector(p-K_neighbor:p+K_neighbor),percent*100); % do the median mask which is equvalent to a 50percentile mask
end
rmmissing(percent_mask)
not_normalize = percent_mask;  % store a copy of not_normalize percent_mask before calculating dff
for i = 1:length(percent_mask)
    if i+1 <= length(percent_mask)
        percent_mask(i)=percent_mask(i+1)-percent_mask(i);
    else
        percent_mask(i)=0;
    end
end

percent_mask = percent_mask/median(raw_vector); % divided by the median of the entire raw data
figure
plot(percent_mask)
hold on
xlabel('Time Series',FontSize = 18)
ylabel(['DFF of Neuron ',num2str(n)],FontSize = 18)

movavgvector = ones(100,1); % number of data point used for each bin 
num = (1/length(movavgvector))*movavgvector; % the moving average filter
den = 1;
y = filter(num,den,percent_mask); % apply a moving average filter to create a baseline 
 
plot(y,'r')
xlabel('Time Steps')
ylabel('DFF and Baseline')
legend('DFF','MA Baseline')
%% Moving Average for Masked Data(not_normalize)
figure
plot(not_normalize)
hold on 
movavgvector = ones(300,1); % number of data point used for each bin 
num = (1/length(movavgvector))*movavgvector; % the moving average filter
den = 1;
y = filter(num,den,not_normalize); % apply a moving average filter to create a baseline 
 
plot(y,'r')
xlabel('Time Steps')
ylabel('DFF and Baseline')
legend('DFF','MA Baseline')
%% PCA(loading plot of all cells)  and K Mean
pcacoeff = pca(dffmatrix); % get the pca coefficient
[~,score,latent] = pca(dffmatrix); % scores are % variance explained by pca
new_pcacoeff = pcacoeff(:,1:2);
%scatter(new_pcacoeff(:,1),new_pcacoeff(:,2),45,"filled") % plot the pca1,2 coeff for all 21 cells
rng(10)  % For reproducibility
% Do PCA
[idx,C] = kmeans(new_pcacoeff,3);
X = new_pcacoeff;

figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',16)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',16)
plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',16)
plot(C(:,1),C(:,2),'kx','MarkerSize',10,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Cluster 3','Center')

title('PCA 1 and PCA 2 Loading for all 21 Neurons',Fontsize = 14)
xlabel(['PCA 1 ', num2str(varexplain1),' %'],Fontsize = 14)
ylabel(['PCA 2 ', num2str(varexplain2),' %'],Fontsize = 14)

% So The Loading shows that we may have two groups of cells clustering
% along the pca 1 direction. See raw data, cells that are close to each
% other have similar feature? 
%% PCA Score Classification
range = 3; % How Many PCA from PCA1 do you want in include?
knumber = 2; % Select Number of clusters.
rng(10)
[idx,C] = kmeans(score(:,1:range),knumber); 
X = score(:,1:2);
varexplain1 = latent(1)/sum(latent)*100; % variance explain bu pca1
varexplain2 = latent(2)/sum(latent)*100;

figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',16)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',16)
%for 3 cluster
%plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',16)
plot(C(:,1),C(:,2),'kx','MarkerSize',10,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Center')
%for 3 cluster
%legend('Cluster 1','Cluster 2','Cluster 3','Center')


title('PCA 1 and PCA 2 Clustering for all 21 Neurons')
xlabel(['PCA 1 ', num2str(varexplain1),' %'],Fontsize = 14)
ylabel(['PCA 2 ', num2str(varexplain2),' %'],Fontsize = 14)
%raw_empty_mask = conv(empty_mask,raw_vector);

%% If this data doesn't look reasonable with K-Means, So we will try GMM(Gaussian Mixature Model)
rng('default')% fix random state for reproducibility
figure
X = dffmatrix(:,[1 5]);

gm = fitgmdist(X,2);
scatter(X(:,1),X(:,2),20,'.') % Scatter plot with points of size 10
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
fcontour(gmPDF,[-1 1])


title('GMM of PCA 1 and PCA 2');
xlabel(['PCA 1 ', num2str(varexplain1),' %'],Fontsize = 14)
ylabel(['PCA 2 ', num2str(varexplain2),' %'],Fontsize = 14)

%calculate how many data there are in a circle
 center_of_circle=[0 0]; 
 radius=6;
 x0 = score(:,1) - center_of_circle(1);
 y0 = score(:,2) - center_of_circle(2);
 disp 'Number of data in this circle'
 sum(x0.^2+y0.^2-radius^2<0)  % gives you the number of elements in the circl


%% Eigenvalue for PCA
figure
bar(latent)
%% Converting denoised signals to Spike Train

figure

cellnumber = 9;
CaDecreaseRate = 0.1;
x = 0:0.00001:2;
y = sin(6*t); % a sin kernel with 

%rawhighpass = highpass(dffmatrix(100:1000,1),499,1e3);
rawnoisereduced = medfilt1(dffmatrix(:,cellnumber)); % apply a median filter to reduce noise
rawhighpass = highpass(rawnoisereduced,499,1e3);

transformed = conv(rawnoisereduced,y,'same'); % conv
threshold = mad(transformed)/0.6754; % set the threshold for spike identification
subplot(211),plot(transformed)
title(['Votage Activity and Threshold Cell Number: ',num2str(cellnumber)],Fontsize = 14)
hold on 
yline(threshold,'r','LineWidth',2)


% Extract Spike Train and plot the EVENT VS TIME(ONLY RUN THIS CODE ONCE!!!!!)
% set those < threshold = 0 ,and > threshold = 1
for i = 1:length(transformed)
  if transformed(i) <=threshold
  transformed(i)=0;
  else
  transformed(i)=1;
  end
end


subplot(212),bar(transformed,5,'k')
title(['Spike Count in Time(msc) Cell Number: ',num2str(cellnumber)],Fontsize = 14)
