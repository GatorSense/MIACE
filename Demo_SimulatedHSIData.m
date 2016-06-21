clear all; clc; 
%%
%Demo Script: Run MI-SMF and MI-ACE on Simulated Data
% This Demo:
% 1) Generates Simulated Hyperspectral Data using
%   /FUMI/gen_synthetic_data_code and spectral signatures from the ASTER
%   Spectral Library
% 2) Runs MI-ACE and MI-SMF on Training Data
% 3) Scores MI-ACE and MI-SMF on Test Data
%
%
%%
%Add Paths
addpath('FUMI/gen_synthetic_data_code')
addpath('FUMI/synthetic_data')
%%
%Set Parameters
NumReps = 3;

%Data Set Generation Parameters
load('E_truth.mat');
denom = sum(E_truth.*E_truth);
E_truth = E_truth./repmat(denom, [size(E_truth,1), 1]);
E_t = E_truth(:,1);
E_minus = E_truth(:,2:end);

%MI ACE Parameters
parameters.globalBackgroundFlag = 0;
parameters.softmaxFlag = 0;
parameters.posLabel = 1;
parameters.negLabel = 0;
parameters.maxIter = 100;
parameters.initType = 1; 

%Training Data Generation Parameters
NumTotalPoints = 500; 
TotalNumberBags = 50; 
settings(1).num_pbags = round(TotalNumberBags/2);
settings(1).num_nbags = round(TotalNumberBags/2);
settings(1).num_points = round(NumTotalPoints/TotalNumberBags);
settings(1).n_tar = 2;
settings(1).N_b = 1;
settings(1).Pt_mean = [.25];
settings(1).sigma = 3;
settings(1).expect_SdB = 40;

%Testing Data Generation Parameters
test_settings(1).num_pbags = 50;
test_settings(1).num_nbags = 50;
test_settings(1).num_points = 500;
test_settings(1).n_tar = 500;
test_settings(1).N_b = 1;
test_settings(1).Pt_mean = [.15];
test_settings(1).sigma = 10;
test_settings(1).expect_SdB = 30;

%%
%Run Methods and Save Results
results = {}; 
table = []; 

for iter = 1:NumReps
    
    disp(['Run ', num2str(iter), ' of ', num2str(NumReps)]);
   
    [X_test,~,~,labels_point_test,~]=gen_multi_tar_mixed_data(E_t,E_minus,test_settings(1).num_pbags,test_settings(1).num_nbags,test_settings(1).num_points,test_settings(1).n_tar,test_settings(1).N_b,test_settings(1).Pt_mean,test_settings(1).sigma,test_settings(1).expect_SdB);
    labels_point_test=reshape(labels_point_test',1,(test_settings(1).num_pbags+test_settings(1).num_nbags)*test_settings(1).num_points);
    
    for n_settings = 1:length(settings);
        num_pbags = settings(n_settings).num_pbags;
        num_nbags = settings(n_settings).num_nbags;
        num_points = settings(n_settings).num_points;
        n_tar = settings(n_settings).n_tar;
        N_b = settings(n_settings).N_b;
        Pt_mean = settings(n_settings).Pt_mean;
        sigma = settings(n_settings).sigma;
        expect_SdB = settings(n_settings).expect_SdB;
        
        %generate training and testing data
        [X_train,~,labels_bag,~,bag_number]=gen_multi_tar_mixed_data(E_t,E_minus,num_pbags,num_nbags,num_points,n_tar,N_b,Pt_mean,sigma,expect_SdB);
        X_train = X_train/max(sqrt(sum(X_train.^2,1)));
        labels_bag=reshape(labels_bag',1,(num_pbags+num_nbags)*num_points);
        for i = 1:num_pbags+num_nbags
            dataBags{i} = X_train(:,bag_number == i)';
            labels(i) = unique(labels_bag(bag_number == i));
        end
    
        %Run SMF init1
        parameters.methodFlag = 0;
        parameters.initType = 1;
        [results{iter,n_settings}.smf.optDict, ~, results{iter,n_settings}.smf.b_mu, results{iter,n_settings}.smf.sig_inv_half] = miTarget(dataBags, labels, parameters);
        [smf_out] = smf_det(X_test,results{iter,n_settings}.smf.optDict',results{iter,n_settings}.smf.b_mu',results{iter,n_settings}.smf.sig_inv_half'*results{iter,n_settings}.smf.sig_inv_half)';
        [results{iter,n_settings}.smf.xx,results{iter,n_settings}.smf.yy,~,results{iter,n_settings}.smf.auc] = perfcurve(labels_point_test,smf_out,1);

        %Run ACE init1
        parameters.methodFlag = 1;
        parameters.initType = 1;
        [results{iter,n_settings}.ace.optDict, ~, results{iter,n_settings}.ace.b_mu, results{iter,n_settings}.ace.sig_inv_half] = miTarget(dataBags, labels, parameters);
        [ace_out] = ace_det(X_test,results{iter,n_settings}.ace.optDict',results{iter,n_settings}.ace.b_mu',results{iter,n_settings}.ace.sig_inv_half'*results{iter,n_settings}.ace.sig_inv_half)';
        [results{iter,n_settings}.ace.xx,results{iter,n_settings}.ace.yy,~,results{iter,n_settings}.ace.auc] = perfcurve(labels_point_test,ace_out,1);
        
        %Construct Table
        method_list = fieldnames(results{1});
        for j = 1:numel(fieldnames(results{1}))
            value = results{iter,n_settings}.(method_list{j});
            table(j,iter,n_settings) = value.auc;
        end
        
    end
end

%%
%Cler Unnecessary files and plot results
clear ace_out bag_number dataBags denom E_minus E_t expect_SdB i iter j labels labels_bag labels_point_test method_list N_b n_settings n_tar num_nbags num_pbags num_points NumTotalPoints Pt_mean sigma smf_out TotalNumberBags value X_test X_train

labels1 = {};
labels2 = {};
figure(101); clf; 
for i = 1:NumReps
    subplot(1,2,1); plot(results{i}.smf.optDict); labels1{end+1} = ['SMF Target Concept - Run', num2str(i)]; hold on; plot(results{i}.ace.optDict); labels1{end+1} = ['ACE Target Concept - Run', num2str(i)];
    subplot(1,2,2); plot(results{i}.smf.xx, results{i}.smf.yy ); labels2{end+1} = ['SMF ROC - Run', num2str(i)]; hold on; plot(results{i}.ace.xx, results{i}.ace.yy); labels2{end+1} = ['ACE ROC - Run', num2str(i)];
end
subplot(1,2,1); legend(labels1); axis([0 212 -0.2 0.2 ]); xlabel('Band Number'); 
subplot(1,2,2); legend(labels2);axis([0 1 0 1]);  xlabel('Probability of False Alarm'); ylabel('Probability of Detection'); 