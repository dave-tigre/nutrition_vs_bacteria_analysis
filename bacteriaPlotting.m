function [ bac_period ] = bacteriaPlotting( bac_stool_dat, name )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

bac_period = [0 0 0 0];
fs = 1; % frequency of sampling was once every 24 hours (assumption)
bac_firstIndex = 2:70; % data obtained in first collection half
bac_secondIndex = 118:191; % data obtained in second collection half

t = (0:length(bac_stool_dat(bac_firstIndex,1)) - 1)/fs;
% plotting
figure()
subplot(2,3,1)
bac_norm_firstIndex = z_normalization( bac_stool_dat(bac_firstIndex,1) );
plot(bac_stool_dat(bac_firstIndex,2),bac_norm_firstIndex)
% axis([ 0 225 min(donorA_stool_calcium(:,2)) max(donorA_stool_calcium(:,2))])
title1 = ['Donor A Stool ', name,  ' Days 1-70'];
title(title1)
xlabel('Collection Days')
ylabel('Count')

subplot(2,3,2)
[autocor,lags] = xcorr(bac_norm_firstIndex,'coeff');
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
bac_period(1,1) = short;
bac_period(1,2) = long;
hold on
pks = plot(lags(lcsh)/fs,pksh,'or', ...
    lags(lclg)/fs,pklg+0.05,'vk');
hold off
legend(pks,[repmat('Period',[2 1])])
axis tight



[bac_pxx1,f] = periodogram(bac_norm_firstIndex,[],[],fs);
subplot(2,3,3)
plot(f,bac_pxx1)
title3 = ['Donor A Stool ', name,  ' Days 1-70 Periodogram'];
title(title3)
xlabel('Frequency (cycles/day)')
ylabel('Magnitude')

subplot(2,3,4)
bac_norm_secondIndex = z_normalization( bac_stool_dat(bac_secondIndex,1) );
plot(bac_stool_dat(bac_secondIndex,2),bac_norm_secondIndex)
% axis([ 0 225 min(donorA_stool_calcium(:,2)) max(donorA_stool_calcium(:,2))])
title4 = ['Donor A Stool ', name,  ' Days 123-202'];
title(title4)
xlabel('Collection Days')
ylabel('Count')

subplot(2,3,5)
[autocor,lags] = xcorr(bac_norm_secondIndex,'coeff');
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
bac_period(1,3) = short;
bac_period(1,4) = long;
hold on
pks = plot(lags(lcsh)/fs,pksh,'or', ...
    lags(lclg)/fs,pklg+0.05,'vk');
hold off
legend(pks,[repmat('Period',[2 1])])
% legend(pks,[repmat('Period',[1 1])])
axis tight

[bac_pxx2,f] = periodogram(bac_norm_secondIndex,[],[],fs);
subplot(2,3,6)
plot(f,bac_pxx2)
title6 = ['Donor A Stool ', name,  ' Days 123-202 Periodogram'];
title(title6)
xlabel('Frequency (Cycles/Day)')
ylabel('Magnitude')


end

