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
% File:   ensemblesTestingMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on May 19, 2018.
%
%
% Function description:
% Computes MSE for the testing dataset given the ensemble metadata

%%%% NO NEED TO MODIFY THIS FILE

solveFor = "KM"

% minimizes cross-validation MSE, log-likelihood, AIC, and BIC, respectively
naiveBaggingEnsembles = [4, 1000, 983, 535];
% minimizes cross-validation MSE, log-likelihood, and AIC/BIC, respectively
randomBaggingEnsembles = [2556, 349, 4044];
% minimizes cross-validation MSE, log-likelihood/AIC, BIC, and objective function, respectively
greedyBaggingEnsembles = [ 3,  5,  4, 12; % ensembles number
                          29, 22, 23,  1]; % lambda number

% minimizes cross-validation MSE, log-likelihood, AIC/BIC, and objective function, respectively
naiveOptimumWeightsEnsembles = [900, 999, 997, 951; % ensembles number
                                 31,   1,  19,  31]; % lambda number
% minimizes cross-validation MSE, log-likelihood/AIC, BIC, and objective function, respectively
randomOptimumWeightsEnsembles = [390, 756, 756, 892; % ensembles number
                                  31,   1,  20,  31]; % lambda number
% minimizes cross-validation MSE, and log-likelihood/AIC/BIC/objective function, respectively
greedyOptimumWeightsEnsembles = [9935, 9344; % ensembles number
                                   26,    1]; % lambda number
% Results LL?
% algorithm             MSE          avg_std  MSE + avg_std^2  criterion    n sets
% naiveBagging:         2.3518 oC^2; 0.1956 oC; 2.3901 oC^2; for minMSE (   4 sets)
% naiveBagging:         2.6917 oC^2; 0.3889 oC; 2.8429 oC^2; for maxLL  (1000 sets)
% naiveBagging:         2.6890 oC^2; 0.3866 oC; 2.8385 oC^2; for minAIC ( 983 sets)
% naiveBagging:         2.5984 oC^2; 0.3235 oC; 2.7031 oC^2; for minBIC ( 535 sets)
% 
% randomBagging:        5.4482 oC^2; 1.6882 oC; 8.2982 oC^2; for minMSE (   4 sets)
% randomBagging:        4.8162 oC^2; 1.3873 oC; 7.4862 oC^2; for maxLL  (  45 sets)
% randomBagging:        4.9327 oC^2; 1.4097 oC; 6.9200 oC^2; for minAIC (   7 sets)
% randomBagging:        4.9327 oC^2; 1.4097 oC; 6.9200 oC^2; for minBIC (   7 sets)
% 
% greedyBagging:        2.2294 oC^2; 0.8625 oC; 2.9733 oC^2; for minMSE (   3 sets)
% greedyBagging:        2.3190 oC^2; 1.3487 oC; 4.1380 oC^2; for maxLL  (   5 sets)
% greedyBagging:        2.3190 oC^2; 1.3487 oC; 4.1380 oC^2; for minAIC (   5 sets)
% greedyBagging:        2.3346 oC^2; 1.4108 oC; 4.3250 oC^2; for minBIC (   4 sets)
% greedyBagging:        2.3716 oC^2; 1.4630 oC; 4.5120 oC^2; for minObj (  12 sets)
% 
% naiveOptimumWeights:  2.2674 oC^2; 0.7126 oC; 2.7752 oC^2; for minMSE (   3 sets)
% naiveOptimumWeights:  2.3295 oC^2; 0.8106 oC; 2.9866 oC^2; for maxLL  (   5 sets)
% naiveOptimumWeights:  2.3261 oC^2; 0.8120 oC; 2.9854 oC^2; for minAIC (   6 sets)
% naiveOptimumWeights:  2.3261 oC^2; 0.8120 oC; 2.9854 oC^2; for minBIC (   6 sets)
% naiveOptimumWeights:  2.2586 oC^2; 0.7777 oC; 2.8634 oC^2; for minObj (   3 sets)
% 
% randomOptimumWeights: 2.3666 oC^2; 0.6666 oC; 2.8110 oC^2; for minMSE (   2 sets)
% randomOptimumWeights: 2.8334 oC^2; 1.5260 oC; 5.1621 oC^2; for maxLL  (   2 sets)
% randomOptimumWeights: 2.8334 oC^2; 1.5260 oC; 5.1621 oC^2; for minAIC (   2 sets)
% randomOptimumWeights: 2.8312 oC^2; 1.5163 oC; 5.1304 oC^2; for minBIC (   2 sets)
% randomOptimumWeights: 2.3726 oC^2; 0.6859 oC; 2.8431 oC^2; for minObj (   2 sets)
% 
% greedyOptimumWeights: 2.1419 oC^2; 1.4316 oC; 4.1914 oC^2; for minMSE (   3 sets)
% greedyOptimumWeights: 2.1717 oC^2; 1.4438 oC; 4.2563 oC^2; for maxLL  (   3 sets)
% greedyOptimumWeights: 2.1717 oC^2; 1.4438 oC; 4.2563 oC^2; for minAIC (   3 sets)
% greedyOptimumWeights: 2.1717 oC^2; 1.4438 oC; 4.2563 oC^2; for minBIC (   3 sets)
% greedyOptimumWeights: 2.1717 oC^2; 1.4438 oC; 4.2563 oC^2; for minObj (   3 sets)


% Ensemble Testing Time: 17.7312s; MSE (std): Total = 2.2674 oC^2 (0.7126 oC); Ts = 2.0846 oC^2 (0.5457 oC); Th = 2.4502 oC^2 (0.8794 oC)                                 
% Ensemble Testing Time: 29.1235s; MSE (std): Total = 2.3295 oC^2 (0.8106 oC); Ts = 2.0643 oC^2 (0.7053 oC); Th = 2.5947 oC^2 (0.9158 oC)                                 
% Ensemble Testing Time: 53.1432s; MSE (std): Total = 2.3261 oC^2 (0.8120 oC); Ts = 2.0630 oC^2 (0.7034 oC); Th = 2.5893 oC^2 (0.9206 oC)                                 
% Ensemble Testing Time: 20.2773s; MSE (std): Total = 2.2586 oC^2 (0.7777 oC); Ts = 2.0726 oC^2 (0.6003 oC); Th = 2.4447 oC^2 (0.9550 oC)                                   
% For each criterion, select the algorithm with the lowest variance:
% minMSE: 
% maxLL:  
% minAIC: 
% minBIC: 
% minObj: 

% For each algorithm, these are the ensembles with the lowest variance:
% naiveBagging:         2.3518 oC^2 (0.1956 oC; 2.3901 oC^2) for minMSE (  4 sets)
% randomBagging:        4.8162 oC^2 (1.3873 oC; 6.7408 oC^2) for maxLL  (  7 sets)
% greedyBagging:        2.2294 oC^2 (0.8625 oC; 2.9733 oC^2) for minMSE (  3 sets)
% naiveOptimumWeights:  2.2674 oC^2 (0.7126 oC; 2.7752 oC^2) for minMSE (  3 sets)
% randomOptimumWeights: 2.3666 oC^2 (0.6666 oC; 2.8110 oC^2) for minMSE (  2 sets)
% greedyOptimumWeights: 2.1419 oC^2 (1.4314 oC; 4.1908 oC^2) for minMSE (  3 sets)

% Given the ensemble number and type, get the sets and calculate. The function
% might as well just accept inputs (e.g., air temperature).

% Loading saved data
folder_ensembleData = locateFolderPath("ensembleData");
folder_ensemble = locateFolderPath("ensemble");
addpath(folder_ensemble);

load([folder_ensembleData "/MSEs" solveFor ".mat"])

load([folder_ensembleData "/naiveBagging" solveFor ".mat"])
load([folder_ensembleData "/randomBagging" solveFor ".mat"])
load([folder_ensembleData "/greedyBagging" solveFor ".mat"])

load([folder_ensembleData "/naiveOptimumWeight" solveFor ".mat"])
load([folder_ensembleData "/randomOptimumWeight" solveFor ".mat"])
load([folder_ensembleData "/greedyOptimumWeight" solveFor ".mat"])


%% Running for the cases without ensemble
% KM: #Set: 0; Time: 4.2334 s; MSE: Total = 7.4951 oC^2; Ts = 10.3716 oC^2; Th = 4.6187 oC^2
problemsSetsMSE(   0, 1, 0, "KM");
% KM: #Set: 1; Time: 4.7173 s; MSE: Total = 7.7960 oC^2; Ts = 10.6798 oC^2; Th = 4.9123 oC^2
problemsSetsMSE(   1, 1, 0, "KM");
% KM: #Set: 7651; Time: 4.7618 s; MSE: Total = 2.3707 oC^2; Ts = 2.0312 oC^2; Th = 2.7102 oC^2
problemsSetsMSE(7651, 1, 0, "KM");

%% Running for naive bagging
% Ensemble Testing Total machine time: 23.1490s; MSE (std; NLL): Total = 2.3518 oC^2 (0.1956 oC; 2662.1205); Ts = 2.1042 oC^2 (0.1583 oC; 1729.4777); Th = 2.5994 oC^2 (0.2328 oC; 932.6428)
predictionsNaiveBaggingMinMSE = ensemblePrediction( MSE_trainings(1:naiveBaggingEnsembles(1), 2), ...
                                 ones(1, naiveBaggingEnsembles(1))/naiveBaggingEnsembles(1), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "naiveBaggingMinMSE");
% Testing if we should use parallel
% 1 core:  22.8444, 24.4453, 24.1357
% 2 cores: 20.5190, 18.2285, 29.4341
% 3 cores: 20.2917, 20.6650, 28.9820
% 4 cores: 23.3877, 26.0571, 23.5508
% 5 cores: 24.5454, 25.4215, 24.2313
% Doesn't look like there is any difference. That's probably due to the fact that
% the internal sqp solver is already parallelized. So, using 1 core, 2, or 3 seem
% to result in similar run times. Since I don't intend to use this run time as a 
% beanchmark, I'll run the codes with 1 core so I can continue using my computer.

% Ensemble Testing Total machine time: 7835.9385s; MSE (std; NLL): Total = 2.6917 oC^2 (0.3889 oC; 611.2503); Ts = 2.3917 oC^2 (0.3828 oC; 270.3642); Th = 2.9917 oC^2 (0.3950 oC; 340.8861)
predictionsNaiveBaggingMaxLL = ensemblePrediction( MSE_trainings(1:naiveBaggingEnsembles(2), 2), ...
                                 ones(1, naiveBaggingEnsembles(2))/naiveBaggingEnsembles(2), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "naiveBaggingMaxLL");

% Ensemble Testing Total machine time: 7662.1007s; MSE (std; NLL): Total = 2.6890 oC^2 (0.3866 oC; 619.9096); Ts = 2.3883 oC^2 (0.3815 oC; 272.0750); Th = 2.9896 oC^2 (0.3917 oC; 347.8345)
predictionsNaiveBaggingMinAIC = ensemblePrediction( MSE_trainings(1:naiveBaggingEnsembles(3), 2), ...
                                 ones(1, naiveBaggingEnsembles(3))/naiveBaggingEnsembles(3), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "naiveBaggingMinAIC");

% Ensemble Testing Total machine time: 4473.9090s; MSE (std; NLL): Total = 2.5984 oC^2 (0.3235 oC; 912.3351); Ts = 2.2948 oC^2 (0.3186 oC; 401.6364); Th = 2.9019 oC^2 (0.3284 oC; 510.6986)
predictionsNaiveBaggingMinBIC = ensemblePrediction( MSE_trainings(1:naiveBaggingEnsembles(4), 2), ...
                                 ones(1, naiveBaggingEnsembles(4))/naiveBaggingEnsembles(4), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "naiveBaggingMinBIC");

                                 
%% Running for random bagging

% Ensemble Testing Total machine time: 36.7885s; MSE (std; NLL): Total = 5.4482 oC^2 (1.6882 oC; 53.7548); Ts = 2.0164 oC^2 (1.7633 oC; -1.5117); Th = 8.8800 oC^2 (1.6131 oC; 55.2665)
sets = MSE_trainings(randomGroupsMSE{randomBaggingEnsembles(1), 2}, 3);
predictionsRandomBaggingMinMSE = ensemblePrediction( sets, ...
                                 ones(1, length(sets))/length(sets), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "randomBaggingMinMSE");
                                 
% Ensemble Testing Total machine time: 1277.7352s; MSE (std; NLL): Total = 4.8162 oC^2 (1.3873 oC; 54.5141); Ts = 2.1925 oC^2 (1.3537 oC; -1.2071); Th = 7.4399 oC^2 (1.4209 oC; 55.7212)
sets = MSE_trainings(randomGroupsMSE{randomBaggingEnsembles(2), 2}, 3);
predictionsRandomBaggingMaxLL = ensemblePrediction( sets, ...
                                 ones(1, length(sets))/length(sets), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "randomBaggingMaxLL");
                                 
% Ensemble Testing Total machine time: 107.0658s; MSE (std; NLL): Total = 4.9327 oC^2 (1.4097 oC; 55.6513); Ts = 2.0640 oC^2 (1.4023 oC; -2.6683); Th = 7.8015 oC^2 (1.4171 oC; 58.3197)
% AIC and BIC
sets = MSE_trainings(randomGroupsMSE{randomBaggingEnsembles(3), 2}, 3);
predictionsRandomBaggingMinAIC = ensemblePrediction( sets, ...
                                 ones(1, length(sets))/length(sets), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "randomBaggingMinAIC");




%% Running for greedy bagging
% Ensemble Testing Total machine time: 39.9196s; MSE (std; NLL): Total = 2.2294 oC^2 (0.8625 oC; 35.0390); Ts = 1.8752 oC^2 (0.7748 oC; 16.3886); Th = 2.5837 oC^2 (0.9502 oC; 18.6504)
sets = setsGreedyBagging(1:greedyBaggingEnsembles(1,1), 6, greedyBaggingEnsembles(2,1));
predictionsGreedyBaggingMinMSE = ensemblePrediction( sets, ...
                                 ones(1, length(sets))/length(sets), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "greedyBaggingMinMSE");
                                 
% Ensemble Testing Total machine time: 62.6426s; MSE (std; NLL): Total = 2.3190 oC^2 (1.3487 oC; 1.5171); Ts = 2.2156 oC^2 (1.3638 oC; -0.7786); Th = 2.4224 oC^2 (1.3336 oC; 2.2956)
% LL and AIC
sets = setsGreedyBagging(1:greedyBaggingEnsembles(1,2), 6, greedyBaggingEnsembles(2,2));
predictionsGreedyBaggingMaxLL = ensemblePrediction( sets, ...
                                 ones(1, length(sets))/length(sets), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "greedyBaggingMaxLL");
                                 
% Ensemble Testing Total machine time: 47.5168s; MSE (std; NLL): Total = 2.3346 oC^2 (1.4108 oC; 1.4182); Ts = 2.2499 oC^2 (1.4640 oC; -0.5901); Th = 2.4193 oC^2 (1.3577 oC; 2.0083)
sets = setsGreedyBagging(1:greedyBaggingEnsembles(1,3), 6, greedyBaggingEnsembles(2,3));
predictionsGreedyBaggingMinBIC = ensemblePrediction( sets, ...
                                 ones(1, length(sets))/length(sets), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "greedyBaggingMinBIC");
                                 
% Ensemble Testing Total machine time: 162.6403s; MSE (std; NLL): Total = 2.3716 oC^2 (1.4630 oC; 4.7377); Ts = 2.0769 oC^2 (1.5497 oC; -1.9116); Th = 2.6663 oC^2 (1.3763 oC; 6.6493)
sets = setsGreedyBagging(1:greedyBaggingEnsembles(1,4), 6, greedyBaggingEnsembles(2,4));
predictionsGreedyBaggingMinObj = ensemblePrediction( sets, ...
                                 ones(1, length(sets))/length(sets), 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "greedyBaggingMinBObj");




%% Running for naive optimum weight

% Ensemble Testing Total machine time: 30.6437s; MSE (std; NLL): Total = 2.2674 oC^2 (0.7126 oC; 118.9762); Ts = 2.0846 oC^2 (0.5457 oC; 84.2747); Th = 2.4502 oC^2 (0.8794 oC; 34.7015)
% - 1 to convert from set_ID to #set
sets = setsOptimumWeights{6, naiveOptimumWeightsEnsembles(2,1)}{naiveOptimumWeightsEnsembles(1,1)} - 1;
weights = naiveOptimumWeights{6, naiveOptimumWeightsEnsembles(2,1)}{naiveOptimumWeightsEnsembles(1,1)};
predictionsNaiveOptimumWeightMinMSE = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "naiveOptimumWeightMinMSE");
                                 
% Ensemble Testing Total machine time: 68.6505s; MSE (std; NLL): Total = 2.3295 oC^2 (0.8106 oC; 69.4718); Ts = 2.0643 oC^2 (0.7053 oC; 34.5187); Th = 2.5947 oC^2 (0.9158 oC; 34.9531)
% - 1 to convert from set_ID to #set
sets = setsOptimumWeights{6, naiveOptimumWeightsEnsembles(2,2)}{naiveOptimumWeightsEnsembles(1,2)} - 1;
weights = naiveOptimumWeights{6, naiveOptimumWeightsEnsembles(2,2)}{naiveOptimumWeightsEnsembles(1,2)};
predictionsNaiveOptimumWeightMaxLL = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "naiveOptimumWeightMaxLL");
                                 
% Ensemble Testing Total machine time: 87.9092s; MSE (std; NLL): Total = 2.3261 oC^2 (0.8120 oC; 69.2425); Ts = 2.0630 oC^2 (0.7034 oC; 34.8089); Th = 2.5893 oC^2 (0.9206 oC; 34.4336)
% - 1 to convert from set_ID to #set
% AIC and BIC
sets = setsOptimumWeights{6, naiveOptimumWeightsEnsembles(2,3)}{naiveOptimumWeightsEnsembles(1,3)} - 1;
weights = naiveOptimumWeights{6, naiveOptimumWeightsEnsembles(2,3)}{naiveOptimumWeightsEnsembles(1,3)};
predictionsNaiveOptimumWeightMinAIC = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "naiveOptimumWeightMinAIC");
                                 
% Ensemble Testing Total machine time: 40.5102s; MSE (std; NLL): Total = 2.2586 oC^2 (0.7777 oC; 87.0666); Ts = 2.0726 oC^2 (0.6003 oC; 61.8060); Th = 2.4447 oC^2 (0.9550 oC; 25.2606)
% - 1 to convert from set_ID to #set
sets = setsOptimumWeights{6, naiveOptimumWeightsEnsembles(2,4)}{naiveOptimumWeightsEnsembles(1,4)} - 1;
weights = naiveOptimumWeights{6, naiveOptimumWeightsEnsembles(2,4)}{naiveOptimumWeightsEnsembles(1,4)};
predictionsNaiveOptimumWeightMinObj = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "naiveOptimumWeightMinObj");
                                 
                                 
%% Running for random optimum weight

% Ensemble Testing Total machine time: 14.0095s; MSE (std; NLL): Total = 2.3666 oC^2 (0.6666 oC; 134.4342); Ts = 2.1876 oC^2 (0.5320 oC; 100.2025); Th = 2.5456 oC^2 (0.8013 oC; 34.2317)
sets = MSE_trainings(randomGroupsMSEOptimumWeight{randomOptimumWeightsEnsembles(1,1), 2}, 2);
weights = randomOptimumWeights{6, randomOptimumWeightsEnsembles(2,1)}{randomOptimumWeightsEnsembles(1,1)};
predictionsRandomOptimumWeightMinMSE = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "randomOptimumWeightMinMSE");

% Ensemble Testing Total machine time: 28.8666s; MSE (std; NLL): Total = 2.8334 oC^2 (1.5260 oC; 10.7172); Ts = 2.2914 oC^2 (1.6230 oC; 0.0776); Th = 3.3754 oC^2 (1.4290 oC; 10.6396)
% LL and AIC
sets = MSE_trainings(randomGroupsMSEOptimumWeight{randomOptimumWeightsEnsembles(1,2), 2}, 2);
weights = randomOptimumWeights{6, randomOptimumWeightsEnsembles(2,2)}{randomOptimumWeightsEnsembles(1,2)};
predictionsRandomOptimumWeightMaxLL = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "randomOptimumWeightMaxLL");

% Ensemble Testing Total machine time: 28.8666s; MSE (std; NLL): Total = 2.8312 oC^2 (1.5163 oC; 10.7118); Ts = 2.3262 oC^2 (1.6128 oC; 0.3510); Th = 3.3361 oC^2 (1.4199 oC; 10.3609)
sets = MSE_trainings(randomGroupsMSEOptimumWeight{randomOptimumWeightsEnsembles(1,3), 2}, 2);
weights = randomOptimumWeights{6, randomOptimumWeightsEnsembles(2,3)}{randomOptimumWeightsEnsembles(1,3)};
predictionsRandomOptimumWeightMinBIC = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "randomOptimumWeightMinBIC");

% Ensemble Testing Total machine time: 14.7644s; MSE (std; NLL): Total = 2.3726 oC^2 (0.6859 oC; 148.0366); Ts = 2.2047 oC^2 (0.5334 oC; 112.1717); Th = 2.5404 oC^2 (0.8385 oC; 35.8649)
% - 1 to convert from set_ID to #set
sets = MSE_trainings(randomGroupsMSEOptimumWeight{randomOptimumWeightsEnsembles(1,4), 2}, 2);
weights = randomOptimumWeights{6, randomOptimumWeightsEnsembles(2,4)}{randomOptimumWeightsEnsembles(1,4)};
predictionsRandomOptimumWeightMinObj = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "randomOptimumWeightMinObj");
                                 
                                                        
%% Running for greedy optimum weight

% Ensemble Testing Total machine time: 110.7107s; MSE (std; NLL): Total = 2.1419 oC^2 (1.4314 oC; -3.3276); Ts = 1.8632 oC^2 (1.3423 oC; -4.8227); Th = 2.4206 oC^2 (1.5204 oC; 1.4951)
sets = setsGreedyOptimumWeights{6, greedyOptimumWeightsEnsembles(2,1)}{greedyOptimumWeightsEnsembles(1,1)} - 1;
weights = greedyWeights{6, greedyOptimumWeightsEnsembles(2,1)}{greedyOptimumWeightsEnsembles(1,1)};
predictionsGreedyOptimumWeightMinMSE = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "greedyOptimumWeightMinMSE");

% Ensemble Testing Total machine time: 41.9153s; MSE (std; NLL): Total = 2.1717 oC^2 (1.4438 oC; -1.5901); Ts = 1.9002 oC^2 (1.4477 oC; -4.2866); Th = 2.4432 oC^2 (1.4399 oC; 2.6965
% LL, AIC, BIC, and Obj
sets = setsGreedyOptimumWeights{6, greedyOptimumWeightsEnsembles(2,2)}{greedyOptimumWeightsEnsembles(1,2)} - 1;
weights = greedyWeights{6, greedyOptimumWeightsEnsembles(2,2)}{greedyOptimumWeightsEnsembles(1,2)};
predictionsGreedyOptimumWeightMaxLL = ensemblePrediction( sets, ...
                                 weights, 0, ...
                                 "KM", NA, NA, NA, NA, NA, nameEnsemble = "greedyOptimumWeightMaxLL");

                                 
% running the mean dataset for the best ensemble
heat_power_mean = [  0,   0,   0,   0, ...
                    60,  60,  60,  60, ...
                   100, 100, 100, 100, ...
                   160, 160, 160, 160, ...
                   200, 200, 200, 200];
Ta_pen_mean = [16.377, 15.878, 15.645, 15.664, ...
               16.349, 15.812, 15.634, 15.739, ...
               16.206, 15.937, 15.676, 15.925, ...
               16.220, 16.123, 15.635, 16.018, ...
               16.131, 15.816, 15.509, 16.156];
                   
predictionsBestEnsembleMeanDataset = bestEnsemblePrediction(heat_power_mean, ...
                                       Ta_pen_mean, NA, NA, NA, NA, ...
                                       "bestEnsembleMeanDataset");
                                       
% running for the surface prediciton
% Ta_pen 0 to 30 oC and 0 to 500 W
heat_power_surface_initial = 0:500;
Ta_pen_surface_initial = 10:0.1:30;

Ta_pen_surface =             repmat(Ta_pen_surface_initial, 1, length(heat_power_surface_initial));
heat_power_surface = reshape(repmat(heat_power_surface_initial,  length(Ta_pen_surface_initial), 1), 1, []);


predictionsBestEnsembleSurface_0_500_10_30 = bestEnsemblePrediction(heat_power_surface, ...
                                       Ta_pen_surface, NA, NA, NA, NA, ...
                                       "bestEnsembleSurface_0_500_10_35");
% calculating energy balance
energyBalance("bestEnsembleSurface_0_500_10_30", Ta_pen_surface_initial, heat_power_surface_initial);
                                       
% Ta_brooder 0 to 40 oC and 0 to 1000 W
% have to run this
heat_power_brooder_surface_initial = 0:1000;
Ta_brooder_surface_initial = 0:0.1:40;

Ta_brooder_surface =         repmat(Ta_brooder_surface_initial, 1, length(heat_power_brooder_surface_initial));
heat_power_brooder_surface = reshape(repmat(heat_power_brooder_surface_initial,  length(Ta_brooder_surface_initial), 1), 1, []);


predictionsBestEnsembleSurfaceBrooder_0_1000_0_40 = bestEnsemblePrediction(heat_power_brooder_surface, ...
                                       NA, Ta_brooder_surface, NA, NA, NA, ...
                                       "bestEnsembleSurfaceBrooder_0_1000_0_40");
% calculating energy balance
energyBalance("bestEnsembleSurfaceBrooder_0_1000_0_40", Ta_brooder_surface_initial, heat_power_brooder_surface_initial);
                                       
% Ta_brooder 0 to 40 oC and 0 to 1000 W.
% This models uses Ta_brooder as an input and adds to it how much the supplemental
% heat would increase Ta_brooder. Hence, heat fluxes calculated using this models
% do not reflect heat fluxes values for given Ta_brooder and heat. They do reflect
% predicted heat flux values for values of Ta_brooder after added by the effect of
% supplemental heat. That is, for supplemental heat > 0, Ta_brooder used in the 
% calculations is > Ta_brooder (Ta_brooder_real = Ta_brooder + effect_of_supp_heat).
heat_power_brooder_added_inc_surface_initial = 0:1000;
Ta_brooder_added_inc_surface_initial = 0:0.1:40;

Ta_brooder_added_inc_surface =         repmat(Ta_brooder_added_inc_surface_initial, 1, length(heat_power_brooder_added_inc_surface_initial));
heat_power_brooder_added_inc_surface = reshape(repmat(heat_power_brooder_added_inc_surface_initial,  length(Ta_brooder_added_inc_surface_initial), 1), 1, []);


predictionsBestEnsembleSurfaceBrooderAddedInc_0_1000_0_40 = bestEnsemblePrediction(heat_power_brooder_added_inc_surface, ...
                                       NA, Ta_brooder_added_inc_surface, NA, NA, NA, ...
                                       "bestEnsembleSurfaceBrooderAddedInc_0_1000_0_40", ...
                                       3, 1, 1);
% calculating energy balance
energyBalance("bestEnsembleSurfaceBrooderAddedInc_0_1000_0_40", ...
              Ta_brooder_added_inc_surface_initial, ...
              heat_power_brooder_added_inc_surface_initial);