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
stimsize = 2000; % number of stimuli per condition '
reliableX = [1 1 0 0 0 0]'; % 1 = reliable, 0 = not reliable
successX = [1 0 0 0 2 2]'; % 1 = successful firing, 0 = unsucess
trialsize= 8; % number of trials on each stimulus
multiple_peaksX = [0 0 0 1 2 2]'; % generate multiple peaks on trace or not
typeX = [0 0 0 0 1 2]' ;
% THERE ARE 4 CONDITIONS IN TOTAL: RELIABLE SUCCESS, RELIABLE FAILURE, UNREALIABLE SINGLEPEAKS, UNRELIABLE MULTIPEAKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for h = 1:6
reliable = reliableX(h);
success = successX(h);
multiple_peaks = multiple_peaksX(h);
type = typeX(h);


gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
 
for j = 1:stimsize
name = ['stimulus',num2str(j)];
   if reliable == 1 && success == 0 && multiple_peaks == 0 && type == 0

    index = 0; 
    amp = abs(normrnd(1,2,[trialsize,1]));
    vo = 1; 
    x = linspace(0,26,26);
    mu = normrnd(15,0.1,[trialsize,1]);
    sig = 4*ones(trialsize,1); 
    noise = 10; noisetype = 'sinusoid'; frequency = 1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak
    
   elseif reliable == 1 && success == 1 && multiple_peaks == 0 && type == 0
    index = 0; 
    amp = abs(normrnd(20,5,[trialsize,1]));
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,0.3,[trialsize,1]);
    sig = 4*ones(trialsize,1);
    noise = 10; noisetype = 'sinusoid'; frequency = 1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak
    


  elseif reliable == 0 && success == 0 && multiple_peaks == 0 && type == 0

    index = 0; 
    %amp = normrnd(26,5,8);
    amp = 4*chi2rnd(2,[trialsize,1]); % 2 degree of freedom chi2 distribution 
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,5,[trialsize,1]); % set the mu to be rand around 15
    sig = 4*rand(trialsize,1);
    noise = 50; noisetype = 'sinusoid'; frequency = 0.1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak
    

  elseif reliable == 0 && success == 0 && multiple_peaks == 1 && type == 0
    index = 1 ; % set the index = 1 to initialize the term. 
    %amp = normrnd(26,5,8);
    amp = 4*chi2rnd(2,[trialsize,1]); % 2 degree of freedom chi2 distribution 
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,5,[trialsize,1]); % std of mu = 5
    sig = 3*rand(trialsize,1);
    noise = 50; noisetype = 'sinusoid'; frequency = 0.1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak 
    
  elseif reliable == 0 && success == 2 && multiple_peaks == 2 && type == 1
    index = 0; 
    amp = 25*rand([trialsize,1]);
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,0.3,[trialsize,1]);
    sig = 4*ones(trialsize,1);
    noise = 10; noisetype = 'sinusoid'; frequency = 1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak
    
   elseif reliable == 0 && success == 2 && multiple_peaks == 2 && type == 2
    index = 0; 
    amp = 25*rand([trialsize,1]);
    vo = 0; 
    x = linspace(0,26,26);
    mu = normrnd(15,0.3,[trialsize,1]);
    sig = 4*ones(trialsize,1);
    noise = 10; noisetype = 'sinusoid'; frequency = 1;
    amp2 = 2.6*chi2rnd(2,[trialsize,1]); % amplitude for the second peak
    mu2 = normrnd(20,5,[trialsize,1]); % set the second peak mu to be no around 20
    sig2 = 3*rand(trialsize,1); % standard deviation of the second peak
  
  end

  if reliable == 1 && success == 0 && multiple_peaks == 0 && type == 0
      disp 'reliable_failure'
      simulateData.reliable_failure.(name) = struct('Mean_Normal',num2cell(mu),'Std_Constant',num2cell(sig),'Amplitude_Normal_Absolute',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace',[], ...
      'PeakLoc',[]);      
    

  elseif reliable == 1 && success == 1 && multiple_peaks == 0 && type == 0
      disp 'reliable_success'
      simulateData.reliable_success.(name) = struct('Mean_Normal',num2cell(mu),'Std_Constant',num2cell(sig),'Amplitude_Normal_Absolute',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace',[], ...
      'PeakLoc',[]);
  
  
  elseif reliable == 0 && success == 0 && multiple_peaks == 0 && type == 0
      disp 'unreliable_singlepeak'
      simulateData.unreliable_singlepeak.(name) = struct('Mean_Normal',num2cell(mu),'Std_Random',num2cell(sig),'Amplitude_Chi2',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace', [] ,...
      'PeakLoc',[]);
        
         
  elseif reliable == 0 && success == 0 && multiple_peaks == 1 && type == 0
      disp 'unreliable_multipeak'
      simulateData.unreliable_multipeak.(name) = struct('Mean_Normal',num2cell(mu),'Std_Random',num2cell(sig),'Amplitude_Chi2',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace',[], ...
      'PeakLoc',[]);

  elseif reliable == 0 && success == 2 && multiple_peaks == 2 && type == 1
      disp 'type1'
      simulateData.unreliable_type1.(name) = struct('Mean_Normal',num2cell(mu),'Std_Random',num2cell(sig),'Amplitude_Chi2',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace',[], ...
      'PeakLoc',[]);

  elseif reliable == 0 && success == 2 && multiple_peaks == 2 && type == 2
      disp 'type2'
      simulateData.unreliable_type2.(name) = struct('Mean_Normal',num2cell(mu),'Std_Random',num2cell(sig),'Amplitude_Chi2',num2cell(amp) ...
      ,'NoiseType',{{'noisetype = sinusoid',['NoiseAmp = ', num2str(noise)],['frequency = ', num2str(frequency)]}},'trial_trace',[], ...
      'PeakLoc',[]);

  end

for i = 1:trialsize
  y = (gaus(x,mu(i),sig(i),amp(i),vo) + index*gaus(x,mu2(i),sig2(i),amp2(i),vo) + noise*sin(frequency*rand(1,26)));
  %Plot gaussian
  %subplot(2,1,j);
  %plot(x, y, 'k-', 'LineWidth',1); axis square; 
  %ylim([0,32])
  %xlim([0,26])
  %hold on ;
 
  loc = maxk(ystore(i,:),2,2);
  loc = {find(ystore(i,:) == loc(1)), loc(1);find(ystore(i,:) == loc(2)), loc(2)};

  if reliable == 1 && success == 1 && multiple_peaks == 0 && type == 0
  simulateData.reliable_success.(name)(i).PeakLoc = loc;
  simulateData.reliable_success.(name)(i).trial_trace = y;

  elseif reliable == 1 && success == 0 && multiple_peaks == 0 && type == 0
  simulateData.reliable_failure.(name)(i).PeakLoc = loc;
  simulateData.reliable_failure.(name)(i).trial_trace = y; 

  elseif reliable == 0 && success == 0 && multiple_peaks == 0 && type == 0
  simulateData.unreliable_singlepeak.(name)(i).PeakLoc = loc;
  simulateData.unreliable_singlepeak.(name)(i).trial_trace = y;

  elseif reliable == 0 && success == 0 && multiple_peaks == 1 && type == 0
  simulateData.unreliable_multipeak.(name)(i).PeakLoc = loc;
  simulateData.unreliable_multipeak.(name)(i).trial_trace = y;

  elseif reliable == 0 && success == 2 && multiple_peaks == 2 && type == 1
  simulateData.unreliable_type1.(name)(i).PeakLoc = loc;
  simulateData.unreliable_type1.(name)(i).trial_trace = y;

  elseif reliable == 0 && success == 2 && multiple_peaks == 2 && type == 2
  simulateData.unreliable_type2.(name)(i).PeakLoc = loc;
  simulateData.unreliable_type2.(name)(i).trial_trace = y;

  end

end
end
end



% second for loop

%first for loop

%%
figure
store = simulateData.reliable_success.stimulus2.trial_trace;
imagesc(corr(store'));
colorbar
% SOME THOUGHTS:
figure
histogram(store(3,:))

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
% 7. Deviation within a cluster(disperstion of peaks' location)
% 8. 

