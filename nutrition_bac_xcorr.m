function [ plots ] = nutrition_bac_xcorr(nutrition_name, nutrition_data, bac_name,bac_data )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

stool_firstIndex = 2:70; % data obtained in first collection half
stool_secondIndex = 118:191; % data obtained in second collection half

bac_firstIndex = 2:70; % data obtained in first collection half
bac_secondIndex = 118:191; % data obtained in second collection half



s1_first = z_normalization( nutrition_data(stool_firstIndex,2) );
s2_first = z_normalization( bac_data(bac_firstIndex,1) );

s1_second = z_normalization( nutrition_data(stool_secondIndex,2) );
s2_second = z_normalization( bac_data(bac_secondIndex,1) );

[xcor1,lag1] = xcorr(s1_first,s2_first,21);
[xcor2,lag2] = xcorr(s1_second,s2_second,21);

figure()
subplot(1,2,1)
plot(lag1,xcor1)
grid on
title1 = [nutrition_name, ' vs ', bac_name, ' Correlation Days 1-70'];
title(title1)
ylabel('Correlation')
xlabel('Lag')


subplot(1,2,2)
plot(lag2,xcor2)
grid on
title1 = [nutrition_name, ' vs ', bac_name, ' Correlation Days 123-202'];
title(title1)
ylabel('Correlation')
xlabel('Lag')
plots=0;
end

