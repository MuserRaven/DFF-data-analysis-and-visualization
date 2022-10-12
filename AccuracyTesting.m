function [momx] = AccuracyTesting(score,P)
simdatafiles = dir('C:\Ben Lab\Project1\simulatedDataofDifferrent Noise Type-20221012T000904Z-001\simulatedDataofDifferrent Noise Type\*.mat');
mom = nan(6,1);
accuracy = nan(6,5);
for u = 1:6
    load(simdatafiles(u).name)
%mean of mean 
% Modify(nuding 
%  1 = reliable success
%  2 = reliable failure
%  3 = unreliable
%score = 0.5 ; P = 0.625; % default
re = {'reliable_success', 'reliable_failure', 'unreliable_singlepeak', 'unreliable_multipeak', 'unreliable_type1'};
store = nan(8,26);
for h = 1:5
prediction = nan(1000,1);
for j = 1:1000
for i  = 1:8
store(i,:) = simulateData.(char(re(h))).(['stimulus',num2str(j)])(i).trial_trace;
%plot(1:26, simulateData.unreliable_multipeak.(['stimulus',num2str(j)])(i).trial_trace,'k')
hold on 
end
ChatterjeeCorr(store);
count = max(store,[],2)>17;
count = sum(count)/8;
aver = (count + ChatterjeeCorr(store))/2;
if aver >= score && count >= P
    disp 'reliable_success'
    prediction(j) = 1;
elseif count == 0 
    disp 'reliable_failure'
    prediction(j) = 2;
else
    disp 'unreliable'
    prediction(j) = 3;

end
end

if strcmp(char(re(h)),'reliable_success')
accuracy(u,1) = length(find(prediction == 1))/length(prediction);

elseif strcmp(char(re(h)),'reliable_failure')
accuracy(u,2) = length(find(prediction == 2))/length(prediction);

elseif strcmp(char(re(h)),'unreliable_singlepeak')
accuracy(u,3) = length(find(prediction == 3))/length(prediction);

elseif strcmp(char(re(h)),'unreliable_multipeak')
accuracy(u,4) = length(find(prediction == 3))/length(prediction);

elseif strcmp(char(re(h)),'unreliable_type1')
accuracy(u,5) = length(find(prediction == 3))/length(prediction);


end
end
mom(u) = mean(accuracy(u,:));
end
momx = mean(mom);
accuracy




