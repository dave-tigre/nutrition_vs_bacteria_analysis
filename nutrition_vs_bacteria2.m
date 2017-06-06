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

stool_Fs = 1*7; % frequency of sampling was once every 24 hours (assumption)
stool_firstIndex = 2:70; % data obtained in first collection half
stool_secondIndex = 118:191; % data obtained in second collection half


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

t = (0:length(donorA_stool_calcium(stool_firstIndex,2)) - 1)/stool_Fs;
% plotting
figure()
subplot(2,2,1)
plot(donorA_stool_calcium(stool_firstIndex,1),donorA_stool_calcium(stool_firstIndex,2))
axis([ 0 225 min(donorA_stool_calcium(:,2)) max(donorA_stool_calcium(:,2))])
title('Donor A Stool Calcium First Days')
xlabel('Collection Days')
ylabel('Calcium (s)')

subplot(2,2,2)
[autocor,lags] = xcorr(donorA_stool_calcium(stool_firstIndex,2),'coeff');
plot(lags/stool_Fs,autocor)
xlabel('Lag (days)')
ylabel('Autocorrelation')
[pksh,lcsh] = findpeaks(autocor);
short = mean(diff(lcsh))/stool_Fs;
hold on
pks = plot(lags(lcsh)/stool_Fs,pksh,'or');
hold off
legend(pks,[repmat('Period: ',[1 1]) num2str([short])])



subplot(2,2,3)
plot(donorA_stool_calcium(stool_secondIndex,1),donorA_stool_calcium(stool_secondIndex,2))
axis([ 0 225 min(donorA_stool_calcium(:,2)) max(donorA_stool_calcium(:,2))])
title('Donor A Stool Calcium Second Half')
xlabel('Collection Days')
ylabel('Calcium (s)')

subplot(2,2,4)
[autocor,lags] = xcorr(donorA_stool_calcium(stool_secondIndex,2),'coeff');
plot(lags/stool_Fs,autocor)
xlabel('Lag (days)')
ylabel('Autocorrelation')
[pksh,lcsh] = findpeaks(autocor);
short = mean(diff(lcsh))/stool_Fs;
hold on
pks = plot(lags(lcsh)/stool_Fs,pksh,'or');
hold off
legend(pks,[repmat('Period: ',[1 1]) num2str([short])])

%%
win = donorA_stool_calcium(stool_secondIndex,2);
h = [1/2 1/2];
binomialCoeff = conv(h,h);
% for n = 1:1
%     binomialCoeff = conv(binomialCoeff,h);
% end

figure
days = 0:length(win)-1;
fDelay = (length(binomialCoeff)-1)/2;
binomialMA = filter(binomialCoeff, 1, win);
plot(days,win, ...
     days-fDelay,binomialMA)

 n = length(binomialMA);
y = fft(binomialMA);
Y = fftshift(y);
fshift = (-n/2:n/2-1)*(stool_Fs/n); % zero-centered frequency range
figure
plot(fshift,Y)


% figure
% subplot(2,1,1)
% plot(donorA_stool_calorie(:,1),donorA_stool_calorie(:,2))
% axis([ 0 225 min(donorA_stool_calorie(:,2)) max(donorA_stool_calorie(:,2))])
% title('Donor A Stool Calorie')
% xlabel('Collection Days')
% ylabel('Calorie (s)')
% 
% figure()
% subplot(2,1,1)
% plot(donorA_stool_carb(:,1),donorA_stool_carb(:,2))
% axis([ 0 225 min(donorA_stool_carb(:,2)) max(donorA_stool_carb(:,2))])
% title('Donor A Stool Carb')
% xlabel('Collection Days')
% ylabel('Carb (s)')
% 
% figure
% subplot(2,1,1)
% plot(donorA_stool_cholesterol(:,1),donorA_stool_cholesterol(:,2))
% axis([ 0 225 min(donorA_stool_cholesterol(:,2)) max(donorA_stool_cholesterol(:,2))])
% title('Donor A Stool Cholesterol')
% xlabel('Collection Days')
% ylabel('Cholesterol (s)')
% 
% figure()
% subplot(2,1,1)
% plot(donorA_stool_fat(:,1),donorA_stool_fat(:,2))
% axis([ 0 225 min(donorA_stool_fat(:,2)) max(donorA_stool_fat(:,2))])
% title('Donor A Stool Fat')
% xlabel('Collection Days')
% ylabel('Fat (s)')
% 
% figure()
% subplot(2,1,1)
% plot(donorA_stool_fiber(:,1),donorA_stool_fiber(:,2))
% axis([ 0 225 min(donorA_stool_fiber(:,2)) max(donorA_stool_fiber(:,2))])
% title('Donor A Stool Fiber')
% xlabel('Collection Days')
% ylabel('Fiber (s)')
% 
% figure()
% subplot(2,1,1)
% plot(donorA_stool_protein(:,1),donorA_stool_protein(:,2))
% axis([ 0 225 min(donorA_stool_protein(:,2)) max(donorA_stool_protein(:,2))])
% title('Donor A Stool Protein')
% xlabel('Collection Days')
% ylabel('Protein (s)')
% 
% figure()
% subplot(2,1,1)
% plot(donorA_stool_satfat(:,1),donorA_stool_satfat(:,2))
% axis([ 0 225 min(donorA_stool_satfat(:,2)) max(donorA_stool_satfat(:,2))])
% title('Donor A Stool Saturated Fat')
% xlabel('Collection Days')
% ylabel('Satureated Fat (s)')
% 
% figure()
% subplot(2,1,1)
% plot(donorA_stool_sodium(:,1),donorA_stool_sodium(:,2))
% axis([0 225 min(donorA_stool_sodium(:,2)) max(donorA_stool_sodium(:,2))])
% title('Donor A Stool Sodium')
% xlabel('Collection Days')
% ylabel('Sodium (s)')
% 
% figure()
% subplot(2,1,1)
% plot(donorA_stool_sugar(:,1),donorA_stool_sugar(:,2))
% axis([ 0 225 min(donorA_stool_sugar(:,2)) max(donorA_stool_sugar(:,2))])
% title('Donor A Stool Sugar')
% xlabel('Collection Days')
% ylabel('Sugar (s)')
