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

% plotting
figure()
subplot(2,5,1)
stool_calcium_norm = z_normalization( donorA_stool_calcium(:,2) );
plot(donorA_stool_calcium(:,1),stool_calcium_norm)
% axis([ 0 225 min(donorA_stool_calcium(:,2)) max(stool_calcium_norm)])
title('Donor A Stool Calcium')
xlabel('Collection Days')
ylabel('Calcium (s)')

subplot(2,5,2)
stool_calorie_norm = z_normalization( donorA_stool_calorie(:,2) );
plot(donorA_stool_calorie(:,1),stool_calorie_norm)
axis([ 0 225 min(donorA_stool_calorie(:,2)) max(donorA_stool_calorie(:,2))])
title('Donor A Stool Calorie')
xlabel('Collection Days')
ylabel('Calorie (s)')

subplot(2,5,3)
plot(donorA_stool_carb(:,1),donorA_stool_carb(:,2))
axis([ 0 225 min(donorA_stool_carb(:,2)) max(donorA_stool_carb(:,2))])
title('Donor A Stool Carb')
xlabel('Collection Days')
ylabel('Carb (s)')

subplot(2,5,4)
plot(donorA_stool_cholesterol(:,1),donorA_stool_cholesterol(:,2))
axis([ 0 225 min(donorA_stool_cholesterol(:,2)) max(donorA_stool_cholesterol(:,2))])
title('Donor A Stool Cholesterol')
xlabel('Collection Days')
ylabel('Cholesterol (s)')

subplot(2,5,5)
plot(donorA_stool_fat(:,1),donorA_stool_fat(:,2))
axis([ 0 225 min(donorA_stool_fat(:,2)) max(donorA_stool_fat(:,2))])
title('Donor A Stool Fat')
xlabel('Collection Days')
ylabel('Fat (s)')

subplot(2,5,6)
plot(donorA_stool_fiber(:,1),donorA_stool_fiber(:,2))
axis([ 0 225 min(donorA_stool_fiber(:,2)) max(donorA_stool_fiber(:,2))])
title('Donor A Stool Fiber')
xlabel('Collection Days')
ylabel('Fiber (s)')

subplot(2,5,7)
plot(donorA_stool_protein(:,1),donorA_stool_protein(:,2))
axis([ 0 225 min(donorA_stool_protein(:,2)) max(donorA_stool_protein(:,2))])
title('Donor A Stool Protein')
xlabel('Collection Days')
ylabel('Protein (s)')

subplot(2,5,8)
plot(donorA_stool_satfat(:,1),donorA_stool_satfat(:,2))
axis([ 0 225 min(donorA_stool_satfat(:,2)) max(donorA_stool_satfat(:,2))])
title('Donor A Stool Saturated Fat')
xlabel('Collection Days')
ylabel('Satureated Fat (s)')

subplot(2,5,9)
plot(donorA_stool_sodium(:,1),donorA_stool_sodium(:,2))
axis([0 225 min(donorA_stool_sodium(:,2)) max(donorA_stool_sodium(:,2))])
title('Donor A Stool Sodium')
xlabel('Collection Days')
ylabel('Sodium (s)')

subplot(2,5,10)
plot(donorA_stool_sugar(:,1),donorA_stool_sugar(:,2))
axis([ 0 225 min(donorA_stool_sugar(:,2)) max(donorA_stool_sugar(:,2))])
title('Donor A Stool Sugar')
xlabel('Collection Days')
ylabel('Sugar (s)')

%% Donor A saliva

% Calcium
donorA_saliva_calcium = zeros(saliva_days,2);
donorA_saliva_calorie = zeros(saliva_days,2);
donorA_saliva_carb = zeros(saliva_days,2);
donorA_saliva_cholesterol = zeros(saliva_days,2);
donorA_saliva_fat = zeros(saliva_days,2);
donorA_saliva_fiber = zeros(saliva_days,2);
donorA_saliva_protein = zeros(saliva_days,2);
donorA_saliva_satfat = zeros(saliva_days,2);
donorA_saliva_sodium = zeros(saliva_days,2);
donorA_saliva_sugar = zeros(saliva_days,2);

for i = 1:length(donorA_saliva_loc)
    donorA_saliva_calcium(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_calcium(i,2) = nutrition_calcium(donorA_saliva_loc(i));
    
    donorA_saliva_calorie(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_calorie(i,2) = nutrition_calorie(donorA_saliva_loc(i));
    
    donorA_saliva_carb(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_carb(i,2) = nutrition_carb(donorA_saliva_loc(i));
    
    donorA_saliva_cholesterol(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_cholesterol(i,2) = nutrition_cholesterol(donorA_saliva_loc(i));
    
    donorA_saliva_fat(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_fat(i,2) = nutrition_fat(donorA_saliva_loc(i));
    
    donorA_saliva_fiber(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_fiber(i,2) = nutrition_fiber(donorA_saliva_loc(i));
    
    donorA_saliva_protein(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_protein(i,2) = nutrition_protein(donorA_saliva_loc(i));
    
    donorA_saliva_satfat(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_satfat(i,2) = nutrition_satfat(donorA_saliva_loc(i));
    
    donorA_saliva_sodium(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_sodium(i,2) = nutrition_sodium(donorA_saliva_loc(i));
    
    donorA_saliva_sugar(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_sugar(i,2) = nutrition_sugar(donorA_saliva_loc(i));
end

% plotting
figure()
subplot(2,5,1)
plot(donorA_saliva_calcium(:,1),donorA_saliva_calcium(:,2))
axis([ 0 225 min(donorA_saliva_calcium(:,2)) max(donorA_saliva_calcium(:,2))])
title('Donor A Saliva Calcium')
xlabel('Collection Days')
ylabel('Calcium (s)')

subplot(2,5,2)
plot(donorA_saliva_calorie(:,1),donorA_saliva_calorie(:,2))
axis([ 0 225 min(donorA_saliva_calorie(:,2)) max(donorA_saliva_calorie(:,2))])
title('Donor A Saliva Calorie')
xlabel('Collection Days')
ylabel('Calorie (s)')

subplot(2,5,3)
plot(donorA_saliva_carb(:,1),donorA_saliva_carb(:,2))
axis([ 0 225 min(donorA_saliva_carb(:,2)) max(donorA_saliva_carb(:,2))])
title('Donor A Saliva Carb')
xlabel('Collection Days')
ylabel('Carb (s)')

subplot(2,5,4)
plot(donorA_saliva_cholesterol(:,1),donorA_saliva_cholesterol(:,2))
axis([ 0 225 min(donorA_saliva_cholesterol(:,2)) max(donorA_saliva_cholesterol(:,2))])
title('Donor A Saliva Cholesterol')
xlabel('Collection Days')
ylabel('Cholesterol (s)')

subplot(2,5,5)
plot(donorA_saliva_fat(:,1),donorA_saliva_fat(:,2))
axis([ 0 225 min(donorA_saliva_fat(:,2)) max(donorA_saliva_fat(:,2))])
title('Donor A Saliva Fat')
xlabel('Collection Days')
ylabel('Fat (s)')

subplot(2,5,6)
plot(donorA_saliva_fiber(:,1),donorA_saliva_fiber(:,2))
axis([ 0 225 min(donorA_saliva_fiber(:,2)) max(donorA_saliva_fiber(:,2))])
title('Donor A Saliva Fiber')
xlabel('Collection Days')
ylabel('Fiber (s)')

subplot(2,5,7)
plot(donorA_saliva_protein(:,1),donorA_saliva_protein(:,2))
axis([ 0 225 min(donorA_saliva_protein(:,2)) max(donorA_saliva_protein(:,2))])
title('Donor A Saliva Protein')
xlabel('Collection Days')
ylabel('Protein (s)')

subplot(2,5,8)
plot(donorA_saliva_satfat(:,1),donorA_saliva_satfat(:,2))
axis([ 0 225 min(donorA_saliva_satfat(:,2)) max(donorA_saliva_satfat(:,2))])
title('Donor A Saliva Saturated Fat')
xlabel('Collection Days')
ylabel('Satureated Fat (s)')

subplot(2,5,9)
plot(donorA_saliva_sodium(:,1),donorA_saliva_sodium(:,2))
axis([ 0 225 min(donorA_saliva_sodium(:,2)) max(donorA_saliva_sodium(:,2))])
title('Donor A Saliva Sodium')
xlabel('Collection Days')
ylabel('Sodium (s)')

subplot(2,5,10)
plot(donorA_saliva_sugar(:,1),donorA_saliva_sugar(:,2))
axis([ 0 225 min(donorA_saliva_sugar(:,2)) max(donorA_saliva_sugar(:,2))])
title('Donor A Saliva Sugar')
xlabel('Collection Days')
ylabel('Sugar (s)')


%% Bacteria Abundance
% % Sample ID Matches
otu_c = 1127;
otu_r = 736;
bacteria1 = [];
loc_stool = 1;
loc_saliva = 1;

% [Sample ID, Location in OTU, Corresponding Collection Day]
otu_stool_sample_loc = num2cell(zeros(length(donorA_stool_sampleID),3));
otu_saliva_sample_loc = num2cell(zeros(length(donorA_saliva_sampleID),3));

% find location of the sample ID in the OTU file
for i = 2:otu_r
    for j = 1:length(donorA_stool_sampleID)
       if(not(isempty(cell2mat(strfind(otu_raw(i,1), donorA_stool_sampleID(j,1))))))
        otu_stool_sample_loc(loc_stool,1) = otu_raw(i,1);
        otu_stool_sample_loc(loc_stool,2) = num2cell(i);
        otu_stool_sample_loc(loc_stool,3) = donorA_stool_sampleID(j,2);
        loc_stool = loc_stool+1; 
       end
    end
    for j =1:length(donorA_saliva_sampleID)
       if(not(isempty(cell2mat(strfind(otu_raw(i,1), donorA_saliva_sampleID(j,1))))))
        otu_saliva_sample_loc(loc_saliva,1) = otu_raw(i,1);
        otu_saliva_sample_loc(loc_saliva,2) = num2cell(i);
        otu_saliva_sample_loc(loc_saliva,3) = donorA_saliva_sampleID(j,2);
        loc_saliva = loc_saliva+1;
       end 
    end
end

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
    filename = ['bacteria_data/stool/bacteria_stool_' num2str(k) '.mat'];
    bac_stool = load(filename);
    bac_stool_dat = bac_stool.bac_stool_data.Bacteria;
    %stool_interest_bac(i) = bac_stool_dat;
%     figure()
%     plot(bac_stool_dat(stool_secondIndex,2),bac_stool_dat(stool_secondIndex,1));
    bac_sequence = ['Bacteria Stool Sequence ' num2str(k)];
    title(bac_sequence)
    xlabel('Collection Days')
    ylabel('Bacteria Amount')
end



%% Stool Correlation
Ts = 86400; % sampling period of once a day in seconds
Fs = 1/Ts;
[C1,lag1] = xcorr(bac_stool_dat(stool_secondIndex,1),donorA_stool_calorie(stool_secondIndex,2),'coeff');
figure
plot(lag1/Fs,C1,'k')
ylabel('Amplitude')
grid on
title('Cross-correlation between Bacteria Sample 48 and Stool Calorie ')

[P1,f1] = periodogram(donorA_stool_calorie(stool_secondIndex,2),[],[],Fs,'power');
[P2,f2] = periodogram(bac_stool_dat(stool_secondIndex,1),[],[],Fs,'power');

figure
t = (0:numel(donorA_stool_calorie(stool_secondIndex,2))-1)/Fs;
subplot(2,2,1)
plot(t,donorA_stool_calorie(stool_secondIndex,2),'k')
ylabel('s1')
grid on
title('Time Series')
subplot(2,2,3)
plot(t,bac_stool_dat(stool_secondIndex,1))
ylabel('s2')
grid on
xlabel('Time (secs)')
subplot(2,2,2)
plot(f1,P1,'k')
ylabel('P1')
grid on
axis tight
title('Power Spectrum')
subplot(2,2,4)
plot(f2,P2)
ylabel('P2')
grid on
axis tight
xlabel('Frequency (Hz)')

sig1 = donorA_stool_calorie(stool_secondIndex,2);
sig2 = bac_stool_dat(stool_secondIndex,1);
[Cxy,f] = mscohere(sig1,sig2,[],[],[],Fs);
Pxy     = cpsd(sig1,sig2,[],[],[],Fs);
phase   = -angle(Pxy)/pi*180;
[pks,locs] = findpeaks(Cxy,'MinPeakHeight',0.05);

figure
subplot(2,1,1)
plot(f,Cxy)
title('Coherence Estimate')
grid on

figure
subplot(2,1,1)
plot(f,Cxy)
title('Coherence Estimate')
grid on
hgca = gca;
hgca.XTick = f(locs);
%hgca.YTick = pks;

subplot(2,1,2)
plot(f,phase)
title('Cross-spectrum Phase (deg)')
grid on

hgca.XTick = f(locs); 
% hgca.YTick = round(phase(locs));
xlabel('Frequency (Hz)')

%% Bacteria Saliva

%{
bac_saliva_sample_interest = [1 12 27 33 48];
% for i=2:101
%     bac_saliva_data = getBacteriaData(otu_raw(1,i),otu_saliva_sample_loc,otu_raw,i);
%     filename = ['bacteria_data\saliva\bacteria_saliva_' num2str(i-1) '.mat'];
%     save(filename,'bac_saliva_data');
% end

% Testing with first 10 Bacteria Typse
figure()
for i=1:30
    filename = ['bacteria_data\saliva\bacteria_saliva_' num2str(i) '.mat'];
    bac_saliva = load(filename);
    bac_saliva_dat = bac_saliva.bac_saliva_data.Bacteria;
    figure()
    plot(bac_saliva_dat(:,2),bac_saliva_dat(:,1));
    bac_sequence = ['Bacteria Saliva Sequence ' num2str(i)];
    title(bac_sequence)
    xlabel('Collection Days')
    ylabel('Bacteria Amount')
end

%}

%% Stool Correlation

%% Spectrogram
figure
spectrogram(bac_stool_dat(stool_secondIndex,1),'yaxis') 
figure
spectrogram(donorA_stool_calorie(stool_secondIndex,2),'yaxis') 