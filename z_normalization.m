function [ norm_sig ] = z_normalization( sig )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sig_mean = mean(sig,'omitnan');
sig_std = std(sig,'omitnan');

norm_sig = (sig-sig_mean)/sig_std;

end

