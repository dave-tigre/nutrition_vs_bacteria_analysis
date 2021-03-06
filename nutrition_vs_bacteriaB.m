%% Nutrition Intake VS Bacteria Abundance
% ECE-S436 
% David Tigreros & John Osguthorpe
% 5/7/2017

%% Initilialization
close all; clc;

%% Import Data
[otu_num,otu_txt,otu_raw] = xlsread('OTU.xlsx');
[tax_num,tax_txt,tax_raw] = xlsread('taxatable.xlsx');
[time_num,time_txt,time_raw] = xlsread('TimeSeries_Metadata.xlsx');

%% Time Series Data
collection_days = (time_num(:,1)); % day of which sample was collected
sampleID = time_raw(2:end-1,1:2);
donor = time_raw(2:end-1,4); % donor strings
nutrition_calcium = (time_num(:,4));
nutrition_calorie = (time_num(:,5));
nutrition_carb = (time_num(:,6));
nutrition_cholesterol = (time_num(:,7));
nutrition_fat = (time_num(:,8));
nutrition_fiber = (time_num(:,9));
nutrition_protein = (time_num(:,10));
nutrition_satfat = (time_num(:,11));
nutrition_sodium = (time_num(:,12));
nutrition_sugar = (time_num(:,13));

nutrition_intake = {['calcium'] ['calorie'] ['carb'] ['cholesterol'] ['fat'] ['fiber'] ['protein'] ['satfat'] ['sodium'] ['sugar']};

%% Donor A Sample ID

% Stool
stool_days = 341;
donorA_stool_loc = zeros(stool_days,1);
donorA_stool_sampleID =num2cell(zeros(stool_days,2));
loc = 1;
% find location of donor A stool samples
for i=1:length(donor)
    if(not(isempty(cell2mat(strfind(donor(i), 'DonorA Stool')))))
        donorA_stool_loc(loc) = i;
        donorA_stool_sampleID(loc,:) = sampleID(i,:);
        loc = loc+1;
    end
end

% Saliva
saliva_days = 286;
donorA_saliva_loc = zeros(saliva_days,1);
donorA_saliva_sampleID = num2cell(zeros(saliva_days,2));
loc = 1;
% find location of donor A saliva samples
for i=1:length(donor)
    if(not(isempty(cell2mat(strfind(donor(i), 'DonorA Saliva')))))
        donorA_saliva_loc(loc) = i;
        donorA_saliva_sampleID(loc,:) = sampleID(i,:);
        loc = loc+1;
    end
end
%% Donor A Stool

stool_Fs = 1;
% Initialize Nutrition Intake
donorA_stool_calcium = zeros(stool_days,2); 
donorA_stool_calorie = zeros(stool_days,2);
donorA_stool_carb = zeros(stool_days,2);
donorA_stool_cholesterol = zeros(stool_days,2);
donorA_stool_fat = zeros(stool_days,2);
donorA_stool_fiber = zeros(stool_days,2);
donorA_stool_protein = zeros(stool_days,2);
donorA_stool_satfat = zeros(stool_days,2);
donorA_stool_sodium = zeros(stool_days,2);
donorA_stool_sugar = zeros(stool_days,2);

for i = 1:length(donorA_stool_loc)
    % Calcium
    donorA_stool_calcium(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_calcium(i,2) = nutrition_calcium(donorA_stool_loc(i));
    
    % Calorie
    donorA_stool_calorie(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_calorie(i,2) = nutrition_calorie(donorA_stool_loc(i));
    
    % Carb
    donorA_stool_carb(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_carb(i,2) = nutrition_carb(donorA_stool_loc(i));
    
    % Cholesterol
    donorA_stool_cholesterol(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_cholesterol(i,2) = nutrition_cholesterol(donorA_stool_loc(i));
    
    % Fat
    donorA_stool_fat(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_fat(i,2) = nutrition_fat(donorA_stool_loc(i));
    
    % Nutrition
    donorA_stool_fiber(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_fiber(i,2) = nutrition_fiber(donorA_stool_loc(i));
    
    % Protein
    donorA_stool_protein(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_protein(i,2) = nutrition_protein(donorA_stool_loc(i));
    
    % Saturated Fat
    donorA_stool_satfat(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_satfat(i,2) = nutrition_satfat(donorA_stool_loc(i));
    
    % Sodium
    donorA_stool_sodium(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_sodium(i,2) = nutrition_sodium(donorA_stool_loc(i));
    
    % Sugar
    donorA_stool_sugar(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_sugar(i,2) = nutrition_sugar(donorA_stool_loc(i));
end
t = (0:length(donorA_stool_calcium(118:191,2)) - 1)/stool_Fs;
% plotting
figure()
% subplot(2,5,1)
plot(donorA_stool_calcium(118:191,1),donorA_stool_calcium(118:191,2))
hold on
plot(t,donorA_stool_calcium(118:191,2))
% axis([ 0 225 min(donorA_stool_calcium(:,2)) max(donorA_stool_calcium(:,2))])
title('Donor A Stool Calcium')
xlabel('Collection Days')
ylabel('Calcium (s)')

[autocor,lags] = xcorr(donorA_stool_calcium(118:191,2),3*7*stool_Fs,'coeff');
figure()
plot(lags/stool_Fs,autocor)
xlabel('Lag (days)')
ylabel('Autocorrelation')

hold on
[pksh,lcsh] = findpeaks(autocor);
short = mean(diff(lcsh))/stool_Fs
[pklg,lclg] = findpeaks(autocor, ...
    'MinPeakDistance',ceil(short)*stool_Fs,'MinPeakheight',0.3);
long = mean(diff(lclg))/stool_Fs
hold on
 pks = plot(lags(lcsh)/stool_Fs,pksh,'or')%, ...
%     lags(lclg)/stool_Fs,pklg+0.05,'vk');
 hold off
legend(pks,[repmat('Period: ',[2 1]) num2str([short;long],0)])
title('Donor A Calcium Periodicity')

figure()
[pxx,f] = periodogram(donorA_stool_calcium(118:191,2),[],[],stool_Fs);

plot(f,pxx)
ax = gca;
ax.XLim = [0 10];
xlabel('Frequency (cycles/week)')
ylabel('Magnitude')

% axis([-21 21 -0.4 1.1])

% t = (0:length(donorA_stool_calorie(:,2)) - 1)/stool_Fs;
% figure()
% plot(donorA_stool_calorie(:,1),donorA_stool_calorie(:,2))
% hold on
% plot(t,donorA_stool_calorie(:,2))
% axis([ 0 225 min(donorA_stool_calorie(:,2)) max(donorA_stool_calorie(:,2))])
% title('Donor A Stool Calorie')
% xlabel('Collection Days')
% ylabel('Calorie (s)')


%% Stool Periodicity

%% Bacteria Stool
stool_firstIndex = 2:70;
stool_secondIndex = 118:191;

bac_stoolinterest = [1 12 27 33 48];

% for i=2:101
%     bac_stool_data = getBacteriaData(otu_raw(1,i),otu_stool_sample_loc,otu_raw,i);
%     filename = ['bacteria_data\stool\bacteria_stool_' num2str(i-1) '.mat'];
%     save(filename,'bac_stool_data');
% end

% Testing with first 10 Bacteria Typse

for i=1:length(bac_stoolinterest)
    k = bac_stoolinterest(i);
    filename = ['bacteria_data\stool\bacteria_stool_' num2str(k) '.mat'];
    bac_stool = load(filename);
    bac_stool_dat = bac_stool.bac_stool_data.Bacteria;
    %stool_interest_bac(i) = bac_stool_dat;
%     figure()
%     plot(bac_stool_dat(stool_secondIndex,2),bac_stool_dat(stool_secondIndex,1));
%     bac_sequence = ['Bacteria Stool Sequence ' num2str(k)];
%     title(bac_sequence)
%     xlabel('Collection Days')
%     ylabel('Bacteria Amount')
end


t = (0:length(bac_stool_dat(stool_secondIndex,1)) - 1)/stool_Fs;
figure
plot(t,bac_stool_dat(stool_secondIndex,1))


[autocor,lags] = xcorr(bac_stool_dat(stool_secondIndex,1),3*7*stool_Fs,'coeff');
figure()
plot(lags/stool_Fs,autocor)
xlabel('Lag (days)')
ylabel('Autocorrelation')

hold on
[pksh,lcsh] = findpeaks(autocor);
short = mean(diff(lcsh))/stool_Fs
[pklg,lclg] = findpeaks(autocor, ...
    'MinPeakDistance',ceil(short)*stool_Fs,'MinPeakheight',0.3);
long = mean(diff(lclg))/stool_Fs
hold on
pks = plot(lags(lcsh)/stool_Fs,pksh,'or')%, ...
%     lags(lclg)/stool_Fs,pklg+0.05,'vk');
hold off
legend(pks,[repmat('Period: ',[2 1]) num2str([short;long],0)])
title('Bacteria 48 Periodicity')