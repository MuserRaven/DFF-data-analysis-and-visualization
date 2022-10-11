%% Modify(nuding 
%  1 = reliable success
%  2 = reliable failure
%  3 = unreliable
store = nan(8,26);
prediction = nan(1000,1);
for j = 1:1000
for i  = 1:8
store(i,:) = simulateData.unreliable_type1.(['stimulus',num2str(j)])(i).trial_trace;
%plot(1:26, simulateData.unreliable_multipeak.(['stimulus',num2str(j)])(i).trial_trace,'k')
hold on 
end
ChatterjeeCorr(store);
count = max(store,[],2)>25;
count = sum(count)/8;
aver = (count + ChatterjeeCorr(store))/2;
if aver >= 0.4 && count >= 0.7
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
prediction = rmmissing(prediction);
accuracy = length(find(prediction == 3))/length(prediction)
%%
figure
plot(store')





