%% Modify(nuding 
%  1 = reliable success
%  2 = reliable failure
%  3 = unreliable
re = {'reliable_success', 'reliable_failure', 'unreliable_singlepeak', 'unreliable_multipeak', 'unreliable_type1'};
store = nan(8,26);
accuracy = nan(5,2);
for h = 1:5
prediction = nan(1000,1);
for j = 1:1000
for i  = 1:8
store(i,:) = simulateData.(char(re(h))).(['stimulus',num2str(j)])(i).trial_trace;
%plot(1:26, simulateData.unreliable_multipeak.(['stimulus',num2str(j)])(i).trial_trace,'k')
hold on 
end
ChatterjeeCorr(store);
count = max(store,[],2)>18;
count = sum(count)/8;
aver = (count + ChatterjeeCorr(store))/2;
if aver >= 0.4 && count >= 0.625
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
accuracy(1,1) = length(find(prediction == 1))/length(prediction);

elseif strcmp(char(re(h)),'reliable_failure')
accuracy(2,1) = length(find(prediction == 2))/length(prediction);

elseif strcmp(char(re(h)),'unreliable_singlepeak')
accuracy(3,1) = length(find(prediction == 3))/length(prediction);

elseif strcmp(char(re(h)),'unreliable_multipeak')
accuracy(4,1) = length(find(prediction == 3))/length(prediction);

elseif strcmp(char(re(h)),'unreliable_type1')
accuracy(5,1) = length(find(prediction == 3))/length(prediction);


end
end





