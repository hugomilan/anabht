%
% ANABHT - ANAlytical solver for steady-state BioHeat Transfer problems in 1D
% 
% Copyright (C) 2018 by Cornell University. All Rights Reserved.
% 
% Written by Hugo Fernando Maia Milan.
% 
% Free for educational, research and non-profit purposes.
% Refer to the license file for details.
%
% 
% File:   separateTrainingTesting.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 13, 2018.
%
%
% Function description:
% Separates the dataset positions into training and testing positions.

% obtain training and testing positions
data_piglets % this file generates data matrix
% Columns are:
%  1: Day
%  2: Time_in_degress
%  3: Hour
%  4: Brooder
%  5: Heat_source_power
%  6: Sequence
%  7: Piglet
%  8: Tr
%  9: Ts
% 10: Ta_brooder
% 11: RH
% 12: Tg
% 13: Irradiance
% 14: Th_IR

% we will get, at least, two samples from each combination Hour:Heat_source_power
% for the testing dataset. Will separate in 75/25 = 130/43
lengthTesting = 43;
testingPosition = [];
positions = cell(20,1);
celli = 0;
for h = 4:7
  for p = [0 60 100 160 200]
    celli += 1;
    a = (data(:, 3) == h & data(:, 5) == p & (sum(!isnan(data(:,[10, 12, 8, 13, 9, 14])), 2) == 6) ).*(1:200)';
    % cleaning possible positions
    for ii = 1:200
      if(a(ii))
        positions{celli}(end+1) = a(ii);
      end
    end
    randomPosition = unidrnd(length(positions{celli}));
    testingPosition(end + 1) = positions{celli}(randomPosition);
    positions{celli}(randomPosition) = [];
    randomPosition = unidrnd(length(positions{celli}));
    testingPosition(end + 1) = positions{celli}(randomPosition);
    positions{celli}(randomPosition) = [];
  end
end
while(length(testingPosition) < lengthTesting)
  celli = unidrnd(length(positions));
  randomPosition = unidrnd(length(positions{celli}));
  testingPosition(end + 1) = positions{celli}(randomPosition);
  positions{celli}(randomPosition) = [];
end
% assemblying the training dataset
trainingPosition = (sum(!isnan(data(:,[10, 12, 8, 13, 9, 14])), 2) == 6).*(1:200)';
% removing the data points that are in the testing dataset
trainingPosition(testingPosition) = [];
% removing the points that are NA
ii = 1;
while (ii <= length(trainingPosition))
  if (trainingPosition(ii) == 0)
    trainingPosition(ii) = [];
  else
    ii += 1;
  end
end

% getting cross-validation positions
% cross_validation_n contains all positions but cross_validation_1c, c for complement
cross_validation_1 = zeros(1, length(trainingPosition)*4/5 );
cross_validation_1c = zeros(1, length(trainingPosition)/5);

cross_validation_2 = zeros(1, length(trainingPosition)*4/5 );
cross_validation_2c = zeros(1, length(trainingPosition)/5);

cross_validation_3 = zeros(1, length(trainingPosition)*4/5 );
cross_validation_3c = zeros(1, length(trainingPosition)/5);

cross_validation_4 = zeros(1, length(trainingPosition)*4/5 );
cross_validation_4c = zeros(1, length(trainingPosition)/5);

cross_validation_5 = zeros(1, length(trainingPosition)*4/5 );
cross_validation_5c = zeros(1, length(trainingPosition)/5);


for ii = 1:(length(trainingPosition)/5)
  cross_validation_1c(ii) = trainingPosition( 1 + (ii - 1)*5 );
  cross_validation_2c(ii) = trainingPosition( 2 + (ii - 1)*5 );
  cross_validation_3c(ii) = trainingPosition( 3 + (ii - 1)*5 );
  cross_validation_4c(ii) = trainingPosition( 4 + (ii - 1)*5 );
  cross_validation_5c(ii) = trainingPosition( 5 + (ii - 1)*5 );
  
  cross_validation_1( (1 + 4*(ii - 1) ):(4*ii) ) = [cross_validation_2c(ii) ...
    cross_validation_3c(ii) cross_validation_4c(ii) cross_validation_5c(ii)];
  cross_validation_2( (1 + 4*(ii - 1) ):(4*ii) ) = [cross_validation_1c(ii) ...
    cross_validation_3c(ii) cross_validation_4c(ii) cross_validation_5c(ii)];
  cross_validation_3( (1 + 4*(ii - 1) ):(4*ii) ) = [cross_validation_2c(ii) ...
    cross_validation_1c(ii) cross_validation_4c(ii) cross_validation_5c(ii)];
  cross_validation_4( (1 + 4*(ii - 1) ):(4*ii) ) = [cross_validation_2c(ii) ...
    cross_validation_3c(ii) cross_validation_1c(ii) cross_validation_5c(ii)];
  cross_validation_5( (1 + 4*(ii - 1) ):(4*ii) ) = [cross_validation_2c(ii) ...
    cross_validation_3c(ii) cross_validation_4c(ii) cross_validation_1c(ii)];
end

save("datasets/trainingPosition.mat", "trainingPosition");
save("datasets/testingPosition.mat", "testingPosition");
save("datasets/crossValidationPosition.mat", "cross_validation_1", ...
     "cross_validation_2", "cross_validation_3", ...
     "cross_validation_4", "cross_validation_5", ...
     "cross_validation_1c", ...
     "cross_validation_2c", "cross_validation_3c", ...
     "cross_validation_4c", "cross_validation_5c");
