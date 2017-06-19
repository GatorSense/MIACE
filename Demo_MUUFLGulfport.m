clear all; clc; 
%%
%Demo: Run MI-SMF and MI-ACE on Gulfport Data
% This Demo:
% 1) Runs MI-ACE and MI-SMF on Training MUUFL Gulfport Data
% 2) Scores MI-ACE and MI-SMF on Test Data using Bullwinkle code
%
%
%%
% Set Paths
addpath('MUUFLGulfport/MUUFLGulfportDataCollection/util/');
addpath('MUUFLGulfport/MUUFLGulfportDataCollection/Bullwinkle/');
addpath('MUUFLGulfport/MUUFLGulfportDataCollection/signature_detectors');
data_path_test = 'MUUFLGulfport/MUUFLGulfportDataCollection/muufl_gulfport_campus_w_lidar_1.mat';
data_path_train = 'MUUFLGulfport/MUUFLGulfportDataCollection/muufl_gulfport_campus_3.mat';

%%
% Set Parameters and Target Type
gulfport_allexp_para.remove_band_index=[1:4 69:72];
gulfport_allexp_para.type = 'brown';

parameters.methodFlag = 0;  %Set to 0 for MI-SMF, Set to 1 for MI-ACE
parameters.initType =1 ; %Options: 1, 2, or 3.  InitType 1 is to use best positive instance based on objective function value, type 2 is to select positive instance with smallest cosine similarity with negative instance mean, type 3 is random selection of instance from positive bag
parameters.globalBackgroundFlag = 0;  %Set to 1 to use global mean and covariance, set to 0 to use negative bag mean and covariance
parameters.softmaxFlag = 0; %Set to 0 to use max, set to 1 to use softmax in computation of objective function values of positive bags
parameters.posLabel = 1; %Value used to indicate positive bags, usually 1
parameters.negLabel = 0; %Value used to indicate negative bags, usually 0 or -1
parameters.maxIter = 100; %Maximum number of iterations (rarely used)

positiveBagSize = 5;
scoring_para= { {gulfport_allexp_para.type,[],[],[]} };

%%
%Create Training Data
load(data_path_train)
hsi.Data=double(hsi.Data);
hsi.Data(:,:,gulfport_allexp_para.remove_band_index) = [];
[n_row,n_col,n_dim]=size(hsi.Data);
index = strcmp(hsi.groundTruth.Targets_Type, gulfport_allexp_para.type);
pBags = {};
nBags = {};
target_col=hsi.groundTruth.Targets_colIndices(index);
target_row=hsi.groundTruth.Targets_rowIndices(index);
labels=zeros(n_row,n_col);
shiftAmt = (positiveBagSize - 1)/2;

hsi_data=reshape(shiftdim(hsi.Data,2),n_dim,n_row*n_col);

for j=1:length(target_col)
    labels(target_row(j)-shiftAmt:target_row(j)+shiftAmt,target_col(j)-shiftAmt:target_col(j)+shiftAmt)=1; %generate binary labels
    temp = hsi.Data(target_row(j)-shiftAmt:target_row(j)+shiftAmt, target_col(j)-shiftAmt:target_col(j)+shiftAmt, :);
    pBags{j} = reshape(temp, [positiveBagSize^2, n_dim]);
end
reshape_labels=reshape(labels,1,n_row*n_col);
valid_mask=reshape(hsi.valid_mask,1,n_row*n_col);
valid_labels=reshape_labels(valid_mask);
valid_hsi_data=hsi_data(:,valid_mask);
lVec = valid_labels(:);
nBags{1} = valid_hsi_data(:,lVec == 0)';

dataBags = horzcat(nBags, pBags);
labelsB = horzcat(zeros(1,length(nBags)), ones(1, length(pBags)));

%%
%Run algorithms
load(data_path_test)
hsi.Data(:,:,gulfport_allexp_para.remove_band_index) = [];
hsi_img = double(hsi.Data);

%%%%%MI-SMF
confid_score={};
parameters.initsampleflag = 0;
parameters.sampltpor = 0.1;
[tar_sig,~, B_mu, sig_inv_half] = miTarget(dataBags, labelsB, parameters);
B_mu = zeros(size(B_mu));
inv_B_cov=sig_inv_half'*sig_inv_half;
confid_out= ace_det_local(hsi_img,tar_sig',hsi.valid_mask,B_mu',inv_B_cov);
confid_score{1}=score_hylid_perpixel(hsi,confid_out,scoring_para,'MISMF','det_fig',[],'roc_fig',150);
figure(150);axis([0 0.001 0 1])

%%%%%MI-ACE
parameters.methodFlag = 1;  %Set to 0 for MI-SMF, Set to 1 for MI-ACE
parameters.initsampleflag = 0;
parameters.sampltpor = 0.1;
[tar_sig] = miTarget(dataBags, labelsB, parameters);
confid_out = ace_det_local(hsi_img,tar_sig',hsi.valid_mask,B_mu',inv_B_cov);
confid_score{2} =score_hylid_perpixel(hsi,confid_out,scoring_para,'MIACE','det_fig',[],'roc_fig',250);
figure(250);axis([0 0.001 0 1])











