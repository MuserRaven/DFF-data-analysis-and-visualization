function [corrresponse] = ChatterjeeCorr(data)
% Calculate chatterjee correlation ;

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

pair = sortrows(pair,1); 
[~,p] = sort(pair(:,2),'ascend');
r = 1:length(pair(:,2));
r(p) = r; yrank = r';

sumy = nan(datalength*trialnum-1,1);
for i = 1:datalength*trialnum-1
sumy(i) = abs(yrank(i+1)-yrank(i));
end

%disp 'ChatterjCorrelation = ' 
corrresponse = 1-(((3*sum(sumy))/((datalength*trialnum)^2 -1)));



%
