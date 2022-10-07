
%x is an array of x-values.
%mu is the mean
%sig is the standard deviation 
%amp is the amplitude
%vo is the vertical offset from baseline (positive or negative)
rng(10)  % set a random seed
simulateData = struct; % Build a struct to store parameters informatio
% ONLY RUN THIS SECTION FOR ONCE, OTHERWISE NEW FIELDS WILL OVERWRITE THE OLD !!!

%% Perfectly reliable Case( both mu and sigma are fixed) And not reliable random cases
%
figure
%%%%%%%%%%%%Set the condition for data generator
reliable = 1; % 1 = reliable, 0 = not reliable
success = 1; % 1 = successful firing, 0 = unsucess
trialsize = 20; % number of trials on each stimulus
multiple_peaks = 0; % generate multiple peaks on trace or not
% THERE ARE 4 CONDITIONS IN TOTAL: RELIABLE SUCCESS, RELIABLE FAILURE, UNREALIABLE SINGLEPEAKS, UNRELIABLE MULTIPEAKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
 
for j = 1:10
name = ['stimulus',num2str(j)];
if reliable == 1 

  if success == 0 
    index = 0; normalization = 1;
    amp = abs(normrnd(1,3,[trialsize,1]));
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,0.3,[trialsize,1]);
    sig = 3*ones(trialsize,1);
    noise = 10; noisetype = 'sinusoid'; frequency = 1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak
    

  elseif success == 1 
    index = 0; normalization = 1;
    amp = abs(normrnd(26,3,[trialsize,1]));
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,0.3,[trialsize,1]);
    sig = 3*rand(trialsize,1);
    noise = 10; noisetype = 'sinusoid'; frequency = 1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak
    
  end

elseif reliable == 0 

  if multiple_peaks == 0
    index = 0; normalization = 1;
    %amp = normrnd(26,5,8);
    amp = 2.6*chi2rnd(2,[trialsize,1]); % 2 degree of freedom chi2 distribution 
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,5,[trialsize,1]); % set the mu to be rand around 15
    sig = 3*rand(trialsize,1);
    noise = 50; noisetype = 'sinusoid'; frequency = 0.1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak
    



  elseif multiple_peaks == 1
    index = 1 ; % set the index = 1 to initialize the term. 
    %amp = normrnd(26,5,8);
    normalization = 2; % set the normalization = 2 to initialize the term.
    amp = 2.6*chi2rnd(2,[trialsize,1]); % 2 degree of freedom chi2 distribution 
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,5,[trialsize,1]); % std of mu = 5
    sig = 3*rand(trialsize,1);
    noise = 50; noisetype = 'sinusoid'; frequency = 0.1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak 
    
    
  end   
end

for i = 1:trialsize
  y = (gaus(x,mu(i),sig(i),amp(i),vo) + index*gaus(x,mu2(i),sig2(i),amp2(i),vo) + noise*sin(frequency*rand(1,26)))/normalization;
  % Plot gaussian
  subplot(5,2,j);
  plot(x, y, 'k-', 'LineWidth',1); axis square; axis off;
  ylim([0,32])
  xlim([0,26])
  hold on ;

  if reliable == 1 && success == 0
      simulateData.reliable_failure.(name) = struct('Mean_Normal',num2cell(mu),'Std_Constant',num2cell(sig),'Amplitude_Normal_Absolute',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace',y);
 
  elseif reliable == 1 && success == 1
      simulateData.reliable_success.(name) = struct('Mean_Normal',num2cell(mu),'Std_Constant',num2cell(sig),'Amplitude_Normal_Absolute',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace',y);
  
  elseif reliable == 0 && multiple_peaks == 0 
      simulateData.unreliable_singlepeak.(name) = struct('Mean_Normal',num2cell(mu),'Std_Random',num2cell(sig),'Amplitude_Chi2',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace',y);
  
  elseif reliable == 0 && multiple_peaks == 1
      simulateData.unreliable_multipeak.(name) = struct('Mean_Normal',num2cell(mu),'Std_Random',num2cell(sig),'Amplitude_Chi2',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace',y);
  end

end

end


% SOME THOUGHTS:

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
