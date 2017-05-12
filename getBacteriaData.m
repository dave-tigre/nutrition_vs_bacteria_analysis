function [ bacteria_info ] = getBacteriaData(bacteria_seq, otu_sample_loc, otu_data, bac_seq_col)
% Summary of this function goes here
%   Detailed explanation goes here
day_length = length(otu_sample_loc(:,2));
bacdata_value = zeros(day_length,2);

seq_field = 'Sequence';
seq_value = bacteria_seq;

sampleID_field = 'SampleID';
sampleID_value = otu_sample_loc(:,1);

days_field = 'Day';
days_value = otu_sample_loc(:,3);

bacdata_field = 'Bacteria';

for i=1:day_length
    if(cell2mat(otu_sample_loc(i,2)) > 0)
        bacdata_value(i,1) = cell2mat(otu_data(cell2mat(otu_sample_loc(i,2)),bac_seq_col));
        bacdata_value(i,2) = cell2mat(otu_sample_loc(i,3));
    end
end
bacdata_value = sortrows(bacdata_value,2);
bacteria_info = struct(seq_field,seq_value,sampleID_field,sampleID_value,...
    days_field,days_value, bacdata_field, bacdata_value);

end

