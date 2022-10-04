Spike_triggered_average(dataStruct.dendrite.dff,100,8)
%% Tuning Curve

t = 1:1:3e2;
f = exp(-t./8);%
for h = 1:length(dataStruct.spine)
spinetrimean = nan(33,26);
  for j = 1:33 % average across stimuli

    spinesim = zeros(26,1);
    for i = 1:8 % average across trials
      
      spinesim = spinesim + squeeze(dataStruct.spine(h).responses(j,i,:));
  
      % plot dff of stimulus 1 over all trials 
    end
    spinesim = spinesim/8;
    spinetrimean(j,:) = spinesim;
  end

% spinestrimean store the spine's response to all 33 stimuli in a length 26
% window
figure
plot(mean(spinetrimean(:,1:6),2))
end
%% Extract spine of interest
% Say we are interest in spine 7 perferred stimulus = 14 dffmax = 0.5871
testdata = squeeze(dataStruct.spine(7).responses(14,:,:));
testdata = reshape(testdata,[8*26,1]);
histogram(testdata,208,"EdgeColor","r")
hold on 
testdata = squeeze(dataStruct.spine(7).responses(4,:,:));
testdata = reshape(testdata,[8*26,1]);
histogram(testdata,208,"EdgeColor","b")
%% T

plot(deconv(dataStruct.spine(1).dffStimTrials(1,:),f))