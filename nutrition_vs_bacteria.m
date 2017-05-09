%% Nutrition Intake VS Bacteria Abundance
% ECE-S436 
% David Tigreros & John
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
sample_days = 627;
donorA_loc = zeros(sample_days,1);
loc = 1;
for i=1:length(donor)
    if(not(isempty(cell2mat(strfind(donor(i), 'DonorA')))))
        donorA_loc(loc) = i;
        loc = loc+1;
    end
end

donorA_sampleID = {};
for i = 1:length(donorA_loc)
    donorA_sampleID(i,:) = sampleID(i,:);
end
%% Donor A Stool
% Stool
stool_days = 341;
donorA_stool_loc = zeros(stool_days,1);
loc = 1;

% find location of donor A stool samples
for i=1:length(donor)
    if(not(isempty(cell2mat(strfind(donor(i), 'DonorA Stool')))))
        donorA_stool_loc(loc) = i;
        loc = loc+1;
    end
end

% Calcium
donorA_stool_calcium = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_calcium(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_calcium(i,2) = nutrition_calcium(donorA_stool_loc(i));
end

% Calorie
donorA_stool_calorie = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_calorie(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_calorie(i,2) = nutrition_calorie(donorA_stool_loc(i));
end

% Carb
donorA_stool_carb = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_carb(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_carb(i,2) = nutrition_carb(donorA_stool_loc(i));
end

% Cholesterol
donorA_stool_cholesterol = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_cholesterol(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_cholesterol(i,2) = nutrition_cholesterol(donorA_stool_loc(i));
end

% Fat
donorA_stool_fat = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_fat(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_fat(i,2) = nutrition_fat(donorA_stool_loc(i));
end

% fiber
donorA_stool_fiber = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_fiber(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_fiber(i,2) = nutrition_fiber(donorA_stool_loc(i));
end

% Protein
donorA_stool_protein = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_protein(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_protein(i,2) = nutrition_protein(donorA_stool_loc(i));
end

% Sat Fat
donorA_stool_satfat = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_satfat(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_satfat(i,2) = nutrition_satfat(donorA_stool_loc(i));
end

% Sodium
donorA_stool_sodium = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_sodium(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_sodium(i,2) = nutrition_sodium(donorA_stool_loc(i));
end

% Sugar
donorA_stool_sugar = zeros(stool_days,2);
for i = 1:length(donorA_stool_loc)
    donorA_stool_sugar(i,1) = collection_days(donorA_stool_loc(i));
    donorA_stool_sugar(i,2) = nutrition_sugar(donorA_stool_loc(i));
end

% Plotting
figure()
subplot(2,5,1)
plot(donorA_stool_calcium(:,1),donorA_stool_calcium(:,2))
title('Donor A Stool Calcium')
xlabel('Collection Days')
ylabel('Calcium (s)')

subplot(2,5,2)
plot(donorA_stool_calorie(:,1),donorA_stool_calorie(:,2))
title('Donor A Stool Calorie')
xlabel('Collection Days')
ylabel('Calorie (s)')

subplot(2,5,3)
plot(donorA_stool_carb(:,1),donorA_stool_carb(:,2))
title('Donor A Stool Carb')
xlabel('Collection Days')
ylabel('Carb (s)')

subplot(2,5,4)
plot(donorA_stool_cholesterol(:,1),donorA_stool_cholesterol(:,2))
title('Donor A Stool Cholesterol')
xlabel('Collection Days')
ylabel('Cholesterol (s)')

subplot(2,5,5)
plot(donorA_stool_fat(:,1),donorA_stool_fat(:,2))
title('Donor A Stool Fat')
xlabel('Collection Days')
ylabel('Fat (s)')

subplot(2,5,6)
plot(donorA_stool_fiber(:,1),donorA_stool_fiber(:,2))
title('Donor A Stool Fiber')
xlabel('Collection Days')
ylabel('Fiber (s)')

subplot(2,5,7)
plot(donorA_stool_protein(:,1),donorA_stool_protein(:,2))
title('Donor A Stool Protein')
xlabel('Collection Days')
ylabel('Protein (s)')

subplot(2,5,8)
plot(donorA_stool_satfat(:,1),donorA_stool_satfat(:,2))
title('Donor A Stool Saturated Fat')
xlabel('Collection Days')
ylabel('Satureated Fat (s)')

subplot(2,5,9)
plot(donorA_stool_sodium(:,1),donorA_stool_sodium(:,2))
title('Donor A Stool Sodium')
xlabel('Collection Days')
ylabel('Sodium (s)')

subplot(2,5,10)
plot(donorA_stool_sugar(:,1),donorA_stool_sugar(:,2))
title('Donor A Stool Sugar')
xlabel('Collection Days')
ylabel('Sugar (s)')

%% Donor A saliva
% saliva
saliva_days = 286;
donorA_saliva_loc = zeros(saliva_days,1);
loc = 1;

% find location of donor A saliva samples
for i=1:length(donor)
    if(not(isempty(cell2mat(strfind(donor(i), 'DonorA Saliva')))))
        donorA_saliva_loc(loc) = i;
        loc = loc+1;
    end
end

% Calcium
donorA_saliva_calcium = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_calcium(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_calcium(i,2) = nutrition_calcium(donorA_saliva_loc(i));
end

% Calorie
donorA_saliva_calorie = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_calorie(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_calorie(i,2) = nutrition_calorie(donorA_saliva_loc(i));
end

% Carb
donorA_saliva_carb = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_carb(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_carb(i,2) = nutrition_carb(donorA_saliva_loc(i));
end

% Cholesterol
donorA_saliva_cholesterol = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_cholesterol(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_cholesterol(i,2) = nutrition_cholesterol(donorA_saliva_loc(i));
end

% Fat
donorA_saliva_fat = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_fat(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_fat(i,2) = nutrition_fat(donorA_saliva_loc(i));
end

% fiber
donorA_saliva_fiber = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_fiber(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_fiber(i,2) = nutrition_fiber(donorA_saliva_loc(i));
end

% Protein
donorA_saliva_protein = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_protein(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_protein(i,2) = nutrition_protein(donorA_saliva_loc(i));
end

% Sat Fat
donorA_saliva_satfat = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_satfat(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_satfat(i,2) = nutrition_satfat(donorA_saliva_loc(i));
end

% Sodium
donorA_saliva_sodium = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_sodium(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_sodium(i,2) = nutrition_sodium(donorA_saliva_loc(i));
end

% Sugar
donorA_saliva_sugar = zeros(saliva_days,2);
for i = 1:length(donorA_saliva_loc)
    donorA_saliva_sugar(i,1) = collection_days(donorA_saliva_loc(i));
    donorA_saliva_sugar(i,2) = nutrition_sugar(donorA_saliva_loc(i));
end

% Plotting
figure()
subplot(2,5,1)
plot(donorA_saliva_calcium(:,1),donorA_saliva_calcium(:,2))
axis([ min(donorA_saliva_calcium(:,1)) max(donorA_saliva_calcium(:,1)) min(donorA_saliva_calcium(:,2)) max(donorA_saliva_calcium(:,2))])
title('Donor A Saliva Calcium')
xlabel('Collection Days')
ylabel('Calcium (s)')

subplot(2,5,2)
plot(donorA_saliva_calorie(:,1),donorA_saliva_calorie(:,2))
axis([ min(donorA_saliva_calorie(:,1)) max(donorA_saliva_calorie(:,1)) min(donorA_saliva_calorie(:,2)) max(donorA_saliva_calorie(:,2))])
title('Donor A Saliva Calorie')
xlabel('Collection Days')
ylabel('Calorie (s)')

subplot(2,5,3)
plot(donorA_saliva_carb(:,1),donorA_saliva_carb(:,2))
axis([ min(donorA_saliva_carb(:,1)) max(donorA_saliva_carb(:,1)) min(donorA_saliva_carb(:,2)) max(donorA_saliva_carb(:,2))])
title('Donor A Saliva Carb')
xlabel('Collection Days')
ylabel('Carb (s)')

subplot(2,5,4)
plot(donorA_saliva_cholesterol(:,1),donorA_saliva_cholesterol(:,2))
axis([ min(donorA_saliva_cholesterol(:,1)) max(donorA_saliva_cholesterol(:,1)) min(donorA_saliva_cholesterol(:,2)) max(donorA_saliva_cholesterol(:,2))])
title('Donor A Saliva Cholesterol')
xlabel('Collection Days')
ylabel('Cholesterol (s)')

subplot(2,5,5)
plot(donorA_saliva_fat(:,1),donorA_saliva_fat(:,2))
axis([ min(donorA_saliva_fat(:,1)) max(donorA_saliva_fat(:,1)) min(donorA_saliva_fat(:,2)) max(donorA_saliva_fat(:,2))])
title('Donor A Saliva Fat')
xlabel('Collection Days')
ylabel('Fat (s)')

subplot(2,5,6)
plot(donorA_saliva_fiber(:,1),donorA_saliva_fiber(:,2))
axis([ min(donorA_saliva_fiber(:,1)) max(donorA_saliva_fiber(:,1)) min(donorA_saliva_fiber(:,2)) max(donorA_saliva_fiber(:,2))])
title('Donor A Saliva Fiber')
xlabel('Collection Days')
ylabel('Fiber (s)')

subplot(2,5,7)
plot(donorA_saliva_protein(:,1),donorA_saliva_protein(:,2))
axis([ min(donorA_saliva_protein(:,1)) max(donorA_saliva_protein(:,1)) min(donorA_saliva_protein(:,2)) max(donorA_saliva_protein(:,2))])
title('Donor A Saliva Protein')
xlabel('Collection Days')
ylabel('Protein (s)')

subplot(2,5,8)
plot(donorA_saliva_satfat(:,1),donorA_saliva_satfat(:,2))
axis([ min(donorA_saliva_satfat(:,1)) max(donorA_saliva_satfat(:,1)) min(donorA_saliva_satfat(:,2)) max(donorA_saliva_satfat(:,2))])
title('Donor A Saliva Saturated Fat')
xlabel('Collection Days')
ylabel('Satureated Fat (s)')

subplot(2,5,9)
plot(donorA_saliva_sodium(:,1),donorA_saliva_sodium(:,2))
axis([ min(donorA_saliva_sodium(:,1)) max(donorA_saliva_sodium(:,1)) min(donorA_saliva_sodium(:,2)) max(donorA_saliva_sodium(:,2))])
title('Donor A Saliva Sodium')
xlabel('Collection Days')
ylabel('Sodium (s)')

subplot(2,5,10)
plot(donorA_saliva_sugar(:,1),donorA_saliva_sugar(:,2))
axis([ min(donorA_saliva_sugar(:,1)) max(donorA_saliva_sugar(:,1)) min(donorA_saliva_sugar(:,2)) max(donorA_saliva_sugar(:,2))])
title('Donor A Saliva Sugar')
xlabel('Collection Days')
ylabel('Sugar (s)')


%% Bacteria Abundance
% Sample ID Matches
otu_c = 1127;
otu_r = 736;
bacteria1 = [];
loc = 1;
otu_sample_loc = zeros(length(donorA_sampleID),1);
% cell2mat(strfind(otu_raw(2,1), donorA_sampleID(2)))

% find location of the sample ID in the OTU file
for i = 2:otu_r
    for j = 1:length(donorA_sampleID)
       if(not(isempty(cell2mat(strfind(otu_raw(i,1), donorA_sampleID(j,1))))))
        otu_sample_loc(loc) = i;
        loc = loc+1;
       end 
    end
end

bac1_occurance = zeros(length(otu_sample_loc),1);
for i = 1:length(otu_sample_loc)
    bac1_occurance(i) = 0;
end
