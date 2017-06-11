function [ nutr_period ] = nutritionPlotting( nutrition_data, name )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nutr_period = [0 0 0 0; 0 0 0 0];
fs = 1; % frequency of sampling was once every 24 hours (assumption)
stool_firstIndex = 2:70; % data obtained in first collection half
stool_secondIndex = 118:191; % data obtained in second collection half
yLabel = [name, ' Amount'];
t = (0:length(nutrition_data(stool_firstIndex,2)) - 1)/fs;
% plotting
figure()
subplot(2,3,1)
norm_firstIndex = z_normalization( nutrition_data(stool_firstIndex,2) );
plot(nutrition_data(stool_firstIndex,1),norm_firstIndex)
% axis([ 0 225 min(nutrition_data(:,2)) max(nutrition_data(:,2))])
title1 = ['Donor A Stool ', name,  ' Days 1-70'];
title(title1)
xlabel('Collection Days')
ylabel(yLabel)

subplot(2,3,2)
[autocor,lags] = xcorr(norm_firstIndex,'coeff');
plot(lags/fs,autocor)
title2 = ['Donor A Stool ', name,  ' Days 1-70 Autocorrelation'];
title(title2)
xlabel('Lag (days)')
ylabel('Autocorrelation')
[pksh,lcsh] = findpeaks(autocor);
short = mean(diff(lcsh))/fs;

[pklg,lclg] = findpeaks(autocor, ...
    'MinPeakDistance',(short)*fs,'MinPeakheight',0.1);
long = mean(diff(lclg))/fs;
nutr_period(1,1) = short;
nutr_period(1,2) = long;
hold on
pks = plot(lags(lcsh)/fs,pksh,'or', ...
    lags(lclg)/fs,pklg+0.05,'vk');
hold off
legend(pks,[repmat('Period',[2 1])])
axis tight
% legend(pks,[repmat('Period',[1 1])])

[pxx,f] = periodogram(norm_firstIndex,[],[],fs);
subplot(2,3,3)
plot(f,pxx)
title3 = ['Donor A Stool ', name,  ' Days 1-70 Periodogram'];
title(title3)
xlabel('Frequency (cycles/day)')
ylabel('Magnitude')
axis ([0 0.5 0 16])



subplot(2,3,4)
norm_secondIndex = z_normalization( nutrition_data(stool_secondIndex,2) );
plot(nutrition_data(stool_secondIndex,1),norm_secondIndex)
% axis([ 0 225 min(nutrition_data(:,2)) max(nutrition_data(:,2))])
title4 = ['Donor A Stool ', name,  ' Days 123-202'];
title(title4)
xlabel('Collection Days')
ylabel(yLabel)

subplot(2,3,5)
[autocor,lags] = xcorr(norm_secondIndex,'coeff');
plot(lags/fs,autocor)
title5 = ['Donor A Stool ', name,  ' Days 123-202 Autocorrelation'];
title(title5)
xlabel('Lag (days)')
ylabel('Autocorrelation')
[pksh,lcsh] = findpeaks(autocor);
short = mean(diff(lcsh))/fs;

[pklg,lclg] = findpeaks(autocor, ...
    'MinPeakDistance',(short)*fs,'MinPeakheight',0.1);
long = mean(diff(lclg))/fs;
nutr_period(1,3) = short;
nutr_period(1,4) = long;
hold on
pks = plot(lags(lcsh)/fs,pksh,'or', ...
    lags(lclg)/fs,pklg+0.05,'vk');
hold off
legend(pks,[repmat('Period',[2 1])])
% legend(pks,[repmat('Period',[1 1])])
axis tight

% [pxx,f] = periodogram(stool_calcium_norm_secondIndex,[],nfft,stool_Fs);
% subplot(2,3,6)
% plot(f,pxx,f,6)
% % yPos = 6.41;
% % hold on
% % plot(get(gca,'xlim'), [yPos yPos],'-r'); % Adapts to x limits of current axes
% % hold off
% title6 = ['Donor A Stool ', name,  ' Days 123-202 PSD'];
% title(title6)
% xlabel('Frequency (cycles/day)')
% ylabel('Magnitude')

[pxx,f] = periodogram(norm_secondIndex,[],[],fs);
subplot(2,3,6)
plot(f,pxx)
xlabel('Frequency (Cycles/Day)')
ylabel('Magnitude')
title6 = ['Donor A Stool ', name,  ' Days 123-202 Periodogram'];
title(title6)
axis ([0 0.5 0 16])

plots = 0;
end

