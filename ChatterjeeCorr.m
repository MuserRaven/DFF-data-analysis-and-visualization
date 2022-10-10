function [corrresponse] = ChatterjeeCorr(data)
% Calculate chatterjee correlation ( xlabel is nudge to the right by a random small number between 0.00001:0.0002 ;

datalength = length(data);
trialnum = length(data(:,1));
pair = nan(datalength*trialnum,2);
for n = 1:trialnum
x = 1:datalength;
y = data(n,:);
  if n == 1
  pair(1:datalength,:) = [x',y'];
  else
  pair((n-1)*datalength +1:(n*datalength),:) = [x',y'];
  end
end

%dt = repmat(0.000000001*rand(20,1),8,1);  % initiate smal nudge by a random nuber
%dt = dt(randperm(length(dt)));
pair(:,1) = pair(:,1); %+ %dt;
test = sortrows(pair,1);
[~,ii]=sort(test(:,2),'ascend');
[~,yrank]=sort(ii);

sumy = nan(datalength*trialnum-1,1);
for i = 1:datalength*trialnum-1
sumy(i) = abs(yrank(i+1)-yrank(i));
end
corrresponse = 1-(((3*sum(sumy))/((datalength*trialnum)^2 -1)));



% TRY NUDGING X BY 0.0001 FOR ALL X TO INCREASE THE NUMBER OF RANK 

%
