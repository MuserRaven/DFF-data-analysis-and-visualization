%
%x is an array of x-values.
%mu is the mean
%sig is the standard deviation 
%amp is the (positive or negative)
%vo is the vertical offset from baseline (positive or negative)


%% Perfectly reliable Case( both mu and sigma are fixed) And not reliable random cases
%
figure
%Set the condition for data generator
perfect = 1;
stimulate = struct;
success = 0; 
trialsize = 20;
multiple_peaks = 0;


for j = 1:6
type = 'perfect';


gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

if perfect == 1 

  if success == 0 
    index = 0; normalization = 1;
    amp = normrnd(1,3,trialsize);
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,0.3,trialsize);
    sig = 3*ones(trialsize,1);
    noise = 10; noisetype = 'sinusoid'; frequency = 1;

  elseif success == 1 
    index = 0; normalization = 1;
    amp = normrnd(26,3,trialsize);
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,0.3,trialsize);
    sig = 3*rand(trialsize,1);
    noise = 10; noisetype = 'sinusoid'; frequency = 1;
  end

elseif perfect == 0 

  if multiple_peaks == 0
    index = 0; normalization = 1;
    %amp = normrnd(26,5,8);
    amp = 2.6*chi2rnd(2,[trialsize,1]); % 2 degree of freedom chi2 distribution 
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,5,trialsize); % set the mu to be rand around 15
    sig = 3*rand(trialsize,1);
    noise = 50; noisetype = 'sinusoid'; frequency = 0.1;

  elseif multiple_peaks == 1
    index = 1 ; % set the index = 1 to initialize the term. 
    %amp = normrnd(26,5,8);
    normalization = 2; % set the normalization = 2 to initialize the term.
    amp = 2.6*chi2rnd(2,[trialsize,1]); % 2 degree of freedom chi2 distribution 
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,5,trialsize); % std of mu = 5
    sig = 3*rand(trialsize,1);
    noise = 50; noisetype = 'sinusoid'; frequency = 0.1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(15,5,trialsize); % set the mu to be rand around 15
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak 
  end   
end

for i = 1:trialsize
  y = (gaus(x,mu(i),sig(i),amp(i),vo) + index*gaus(x,mu2(i),sig2(i),amp2(i),vo) + noise*sin(frequency*rand(1,26)))/normalization;
  % Plot gaussian
  subplot(3,2,j);
  plot(x, y, 'k-', 'LineWidth',1); axis square; axis off;
  ylim([0,32])
  xlim([0,26])
  hold on ;



if strcmp(type,'perfect') && strcmp(success,'1')
    struct.perfect.('success')(j) = 1;


    struct.rng= 10 ;
elseif strcmp(type,'perfect') && strcmp(success,'0')
    struct.perfect.('success')(j) = 0;

end

end
end


% Label simulated datasets as 'reliable' or 'unreliable'
% reliable can have two states 1. reliable failure(non-preferred stimulus)
% 2. reliable success(preferred stimulus) 
% Parameters of Interests 1.Correlation Martrix for All Trails in a single
% stimulus, Higher correlation across all trial will means traces are compact
% 2. The Peak Value/Values and where they occur
% 3. The Spatiotemporal distribution of Peak values
% 4. The Distribution/historgram of all trials given stimulus s
% 5. The std of all trials given a stimulus s 
% 6. Klustering properties of peaks values.
% 7. Deviation within a cluster
% 8. 





