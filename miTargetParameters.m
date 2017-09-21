function parameters = miTargetParameters()

% MIACE/MISMF Multiple Instance Adaptive Cosine Estimator/Multiple Instance
%       Spectral Matched Filter Demo
%
% Syntax:  parameters = miTargetParameters()
%
% Outputs:
%   endmembers - double Mat - NxM matrix of M endmembers with N spectral
%       bands
%   P - double Mat - NxM matrix of abundances corresponding to M input
%       pixels and N endmembers
%
% Author: Alina Zare
% University of Florida, Electrical and Computer Engineering
% Email Address: azare@ufl.edu
% Latest Revision: September 21, 2017
% This product is Copyright (c) 2017 University of Florida
% All rights reserved.


parameters.methodFlag = 0;  %Set to 0 for MI-SMF, Set to 1 for MI-ACE
parameters.initType = 1; %Options: 1, 2, or 3.  InitType 1 is to use best positive instance based on objective function value, type 2 is to select positive instance with smallest cosine similarity with negative instance mean, type 3 uses the cluster centers after applying K-means clustering
parameters.globalBackgroundFlag = 0;  %Set to 1 to use global mean and covariance, set to 0 to use negative bag mean and covariance
parameters.softmaxFlag = 0; %Set to 0 to use max, set to 1 to use softmax in computation of objective function values of positive bags (This is generally not used and fixed to 0)
parameters.posLabel = 1; %Value used to indicate positive bags, usually 1
parameters.negLabel = 0; %Value used to indicate negative bags, usually 0 or -1
parameters.maxIter = 1000; %Maximum number of iterations (rarely used)
parameters.samplePor = 1; % If using init1, percentage of positive data points used to initialize (default = 1) 
parameters.initK = 1000; % If using init3, number of clusters used to initialize (default = 1000);
