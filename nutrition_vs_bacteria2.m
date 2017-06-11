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

stool_Fs = 1; % frequency of sampling was once every 24 hours (assumption)
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


%% Nutrition Intake Periodicity

calciumPlots = nutritionPlotting(donorA_stool_calcium,'Calcium');
caloriePlots = nutritionPlotting(donorA_stool_calorie,'Calorie');
carbPlots = nutritionPlotting(donorA_stool_carb,'Carb');
cholesterolPlots = nutritionPlotting(donorA_stool_cholesterol,'Cholesterol');
fatPlots = nutritionPlotting(donorA_stool_fat,'Fat');
fiberPlots = nutritionPlotting(donorA_stool_fiber,'Fiber');
proteinPlots = nutritionPlotting(donorA_stool_protein,'Protein');
satfatPlot = nutritionPlotting(donorA_stool_satfat,'Saturated Fat');
sodiumPlots = nutritionPlotting(donorA_stool_sodium,'Sodium');
sugarPlots = nutritionPlotting(donorA_stool_sugar,'Sugar');

nutrPeriod = [calciumPlots;caloriePlots;carbPlots;cholesterolPlots;fatPlots;fiberPlots;proteinPlots;...
    satfatPlot;sodiumPlots;sugarPlots]

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

bac_stoolinterest = [1 20 33 48 87];
bacPeriods  = [];
eucDis = [];
r = 1;
for i=1:length(bac_stoolinterest)
    k = bac_stoolinterest(i);
    filename = ['bacteria_data/stool/bacteria_stool_' num2str(k) '.mat'];
    bac_stool = load(filename);
    bac_stool_dat = bac_stool.bac_stool_data.Bacteria;
    bac_sequence = ['Bacteria Stool Sequence ' num2str(k)];
    bacName = ['Bacteria ' num2str(k)];
%     bacPeriods(i,:) = bacteriaPlotting(bac_stool_dat,bacName);
%     nutrition_bac_xcorr('Calcium', donorA_stool_calcium, bacName,bac_stool_dat );
%     nutrition_bac_xcorr('Calorie', donorA_stool_calorie, bacName,bac_stool_dat );
%     nutrition_bac_xcorr('Carb', donorA_stool_carb, bacName,bac_stool_dat );
%     nutrition_bac_xcorr('Fat', donorA_stool_fat, bacName,bac_stool_dat );
%     nutrition_bac_xcorr('Fiber', donorA_stool_fiber, bacName,bac_stool_dat );
%     nutrition_bac_xcorr('Sugar', donorA_stool_sugar, bacName,bac_stool_dat );
%    
    figure
eucDis(r,1) = dtw( z_normalization( donorA_stool_calorie(stool_secondIndex,2)),z_normalization( bac_stool_dat(stool_secondIndex,1) ));
    legend('Calorie Intake',bacName)
    xlabel('Collection Days')
    ylabel('Amplitude')
    figure
    eucDis(r,2) = dtw( z_normalization( donorA_stool_calcium(stool_secondIndex,2)),z_normalization( bac_stool_dat(stool_secondIndex,1) ));
     legend('Calcium Intake',bacName)
    xlabel('Collection Days')
    ylabel('Amplitude')
    figure
    eucDis(r,3) = dtw( z_normalization( donorA_stool_carb(stool_secondIndex,2)),z_normalization( bac_stool_dat(stool_secondIndex,1) ));
     legend('Carb Intake',bacName)
    xlabel('Collection Days')
    ylabel('Amplitude')
    figure
    eucDis(r,4) = dtw( z_normalization( donorA_stool_fat(stool_secondIndex,2)),z_normalization( bac_stool_dat(stool_secondIndex,1) ));
     legend('Fat Intake',bacName)
    xlabel('Collection Days')
    ylabel('Amplitude')
    figure
    dtw( z_normalization( donorA_stool_fiber(stool_secondIndex,2)),z_normalization( bac_stool_dat(stool_secondIndex,1) ));
     legend('Fiber Intake',bacName)
    xlabel('Collection Days')
    ylabel('Amplitude')
    figure
    eucDis(r,5) = dtw( z_normalization( donorA_stool_sugar(stool_secondIndex,2)),z_normalization( bac_stool_dat(stool_secondIndex,1) ));
 legend('Sugar Intake',bacName)
    xlabel('Collection Days')
    ylabel('Amplitude')
    
    r = r+1;
end

bacPeriods

%% Cross Correlation

%%
% Fs = 1;
% L = length(bac_norm_secondIndex);
% Y = fft(bac_norm_secondIndex);
% P2 = abs(Y/L);
% P1t = P2(1:L/2+1);
% P1t(2:end-1) = 2*P1t(2:end-1);
% figure
% f = Fs*(0:(L/2))/L;
% plot(f,P1t) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% Fs = 1;
% L = length(z_normalization(donorA_stool_calcium(stool_secondIndex,2)));
% Y2 = fft(z_normalization(donorA_stool_calcium(stool_secondIndex,2)));
% P2 = abs(Y2/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% figure
% f = Fs*(0:(L/2))/L;
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% figure
% corr = P1t.*conj(P1);
% plot(corr)
% 
% 
% figure 
% [corr,lag] = xcorr((z_normalization(donorA_stool_calcium(stool_secondIndex,2))),bac_norm_secondIndex,'coeff');
% plot(corr)
%% Dynamic Time Warping
s1_first = z_normalization( donorA_stool_calorie(stool_secondIndex,2) );
s2_first = z_normalization( bac_stool_dat(stool_secondIndex,1) );
dtw(s1_first,s2_first)

%% more
A = s1_first;
B = s2_first
T = stool_secondIndex;
 X=(ifft(fft(A).*conj(fft(B))));
 figure
 plot(X)
 shiftindex=find(X==max(X));
 shift=T(shiftindex)-T(1)  % This is the time value of the shift.
 
  figure
 plot(A);
 hold on;
 plot(circshift(B,(find(X==max(X)))));
 
 legend('Calorie','Bacteria')
 
 %%
 
 t = 0:0.001:1-0.001;
x = cos(2*pi*100*t);
[xc,lags] = xcorr(x,'coeff');
stem(lags,xc)
stem(lags(length(x):length(x)+50),xc(length(x):length(x)+50));
xlabel('Lags'); ylabel('ACF');