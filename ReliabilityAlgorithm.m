%% isolate spine activity preceeding dendrite spike and extracting timestamp from spine
spine = 5;
meandff = nan(33,1);
for h = 1:33
for n = 1:8
trial1 = squeeze(dataStruct.spine(spine).responses(h,n,:));
%imagesc(trial1)
%ylabel('Stimuli ID')
%xlabel('Presented Time')
%colorbar
%max(trial1)
meandff(n) = mean(trial1,'all');
orientation = squeeze(dataStruct.spine(spine).responses(:,1,1));
end
end
hold on 
orientation = [orientation meandff];
%orientation = sort(orientation,1);
[~,idx] = sort(orientation(:,1)); % sort just the first column
orientation = orientation(idx,:); 
plot(orientation(:,1),orientation(:,2))

%% Demonstrate how algorithm works

rng(10)% set a random seed
figure
sizesim = 500;
dffmax = 6;
%x = normrnd(2,3,[size,1]);
%x = rand([size,1]);
x = random('chi2',1,sizesim,1);

% Here we can assume the dff in non-preferred orientation yields a chi-square 
% or gamma distribution
% And the preferred orientation yields a normal distribution
% This code also works for normal distribution, even it is slightly less
% accurate
x = x(x>=0);
histogram(x,sizesim,"EdgeColor","b")
hold on
y = normrnd(6,1,[sizesim,1]);
y = y(y>=0);
histogram(y,sizesim,"EdgeColor","r")
xlim([0,10])
ylim([0,20])
xlabel('DFF Amplitude')
ylabel('# of Counts')

stdy = 1/(1+exp(-(std(y))));
stdx = 1/(1+exp(-(std(x))));

stdy = std(y);
stdx = std(x);

% But Can we optimize this number?
%stdy and stdx should be both small number 
distance = abs((dffmax-geomean(y))-(dffmax-geomean(x)));
%calculate how far dffmax is from the geomean of x and y
% If we have multiple preferred orientation 
% distance =
% sum(abs(([dffmax1,dffmax2....dffmaxn]-geomean(y1,y2...yn))-(([dffmax1,dffmax2....dffmaxn]-geomean(x1,x2...xn))))/n
   % weight;

   if abs(mean(x)-median(x)) < 0.1 
    disp 'Normal Distribution or Uniform Distribution'
    outlier = find(x(x>dffmax));
    relia = 1/(1+exp(-(distance)/(exp(1.2*stdy)+exp(stdx)))) - (length(outlier)/sizesim);
    % exponential nonlinearity to punish high std 
    % 2 is given to stdy because we want to punish high stdy more than stdx 
    % sigmoid function to normalized the final result
    % calculate the number out outlier in x and punish outlier by
    % outlier/size
    % For multiple distribution, we will take the stdy and stdx for all
    % preffered and non preferred orientation, then compute the mean
    % reliabillity 
   else
    disp 'Non Normal Distribution'
    outlier = find(x(x>dffmax));
    relia = 1/(1+exp(-(distance)/(1.2*exp(stdy)+exp(stdx)))) - (length(outlier)/sizesim);
   end

%%

t = 1:1:3e2;
f = exp(-t./30);
trialtim = nan(12,8,33,26);  % n1 = spine number n2= total number of trial n3 = 
% stimul id, n4 = number of stimuli in each trial
for j = 1:length(dataStruct.spine)

  for n = 1:8

   for i = 1:33  
       
       for k = 1:26
         trialstim(j,n,i,k) = dataStruct.spine(j).responses(i,n,k);
         % this is gonna give out 26 numbers trial n, stimulus i
       end
      
   end

  end

end
trialstim(1,:,:,:)
%plot(orientation(:,1), orientation(:,2))
%%
t = 1:1:3e2;
f = exp(-t./30);
subplot(211)
plot(dffstim)
xline(678)
xline(678*2)
xline(678*3)
xline(678*4)
subplot(212)
plot(dataStruct.spine(1).dff)