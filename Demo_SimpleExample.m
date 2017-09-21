%%
%Demo Script: Run MI-SMF and MI-ACE on Simulated Hyperspectral Data
% This Demo:
% 1) Loads a simulated hyperspectral training and testing set
% 2) Runs MI-ACE and MI-SMF on Training Data
% 3) Scores MI-ACE and MI-SMF on Test Data
%


% MIACE/MISMF Multiple Instance Adaptive Cosine Estimator/Multiple Instance
%       Spectral Matched Filter Demo
%
%
% Author: Alina Zare
% University of Florida, Electrical and Computer Engineering
% Email Address: azare@ufl.edu
% Created: August 2017
% Latest Revision: September 21, 2017
% This product is Copyright (c) 2017 University of Florida
% All rights reserved.

%%
%Set Parameters

%MI ACE Parameters
parameters.globalBackgroundFlag = 0;
parameters.softmaxFlag = 0;
parameters.posLabel = 1;
parameters.negLabel = 0;
parameters.maxIter = 100;
parameters.initType = 1; 

load simpleExampleData; 

%%
%Run Methods and Save Results
    
%Run SMF init1
parameters.methodFlag = 0;
parameters.initType = 1;
parameters.samplePor = 1;
[results.smf.optDict, ~, results.smf.b_mu, results.smf.sig_inv_half] = miTarget(dataBags, labels, parameters);
[smf_out] = smf_det(X_test,results.smf.optDict',results.smf.b_mu',results.smf.sig_inv_half'*results.smf.sig_inv_half)';
[results.smf.xx,results.smf.yy,~,results.smf.auc] = perfcurve(labels_point_test,smf_out,1);

%Run ACE init1
parameters.methodFlag = 1;
parameters.initType = 1;
parameters.samplePor = 1;
[results.ace.optDict, ~, results.ace.b_mu, results.ace.sig_inv_half] = miTarget(dataBags, labels, parameters);
[ace_out] = ace_det(X_test,results.ace.optDict',results.ace.b_mu',results.ace.sig_inv_half'*results.ace.sig_inv_half)';
[results.ace.xx,results.ace.yy,~,results.ace.auc] = perfcurve(labels_point_test,ace_out,1);

%%
%Plot results
labels1 = {};
labels2 = {};
figure(101); clf; 
subplot(1,2,1); plot(results.smf.optDict); labels1{end+1} = ['SMF Target Concept']; hold on; plot(results.ace.optDict); labels1{end+1} = ['ACE Target Concept'];
subplot(1,2,2); plot(results.smf.xx, results.smf.yy ); labels2{end+1} = ['SMF ROC']; hold on; plot(results.ace.xx, results.ace.yy); labels2{end+1} = ['ACE ROC'];
subplot(1,2,1); legend(labels1); axis([0 212 -0.2 0.2 ]); xlabel('Band Number'); 
subplot(1,2,2); legend(labels2);axis([0 1 0 1]);  xlabel('Probability of False Alarm'); ylabel('Probability of Detection'); 

