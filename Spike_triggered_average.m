function [avg, nEvs] = Spike_triggered_average(stim, windowSize)

% magic number for spike
threshold = 5*median(abs(stim)/0.6745);

resp = stim>threshold;
% response
subplot(2,1,1), plot(stim)
line([1, length(stim)], [threshold, threshold], 'LineStyle', '--', ...
    'color', 'k')
ylabel('DFF'), title('Stimulus');
subplot(2,1,2), plot(resp(1:length(stim)), 'r')
ylabel('Spikes'), title('Response')
xlabel('Time')


resp(1:windowSize)=0;

nEvs = sum(resp); % the number of responses
disp 'number of response nEvs', num2str(nEvs)
evIdx = find(resp); % the index of responses

%Save the average in the variable avg, so preallocate it before the for loop.
avg = zeros(1,windowSize);
for w = 1:nEvs
    % Find the indexes of the time window preceding the event
    wIdx = evIdx(w)-windowSize : evIdx(w)-1;
    % Add the dffdendrite from this window to the average
    avg = avg + stim(wIdx);
end

avg=avg./sum(resp);

figure
plot(avg)
title(['Average of ', num2str(nEvs),  'windows'])
xlabel('Time')
ylabel('DFF')


end