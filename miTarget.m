function [optTarget, optObjVal, b_mu, sig_inv_half, init_t] = miTarget(dataBags, labels, parameters)

if(nargin <3)
    parameters.methodFlag = 0;  %Set to 0 for MI-SMF, Set to 1 for MI-ACE
    parameters.globalBackgroundFlag = 0;  %Set to 1 to use global mean and covariance, set to 0 to use negative bag mean and covariance
    parameters.initType = 1; %Options: 1, 2, or 3.  InitType 1 is to use best positive instance based on objective function value, type 2 is to select positive instance with smallest cosine similarity with negative instance mean, type 3 is random selection of instance from positive bag
    parameters.softmaxFlag = 0; %Set to 0 to use max, set to 1 to use softmax in computation of objective function values of positive bags
    parameters.posLabel = 1; %Value used to indicate positive bags, usually 1
    parameters.negLabel = 0; %Value used to indicate negative bags, usually 0 or -1
    parameters.maxIter = 100; %Maximum number of iterations (rarely used)
end

nBags = length(dataBags);
nDim = size(dataBags{1},2);
nPBags = sum(labels == parameters.posLabel);

%Estimate background mean and inv cov
if(parameters.globalBackgroundFlag)    
    data = vertcat(dataBags{:});
    b_mu = mean(data);
    b_cov = cov(data)+eps*eye(size(data,2));
else
    nData = vertcat(dataBags{labels == parameters.negLabel});
    b_mu = mean(nData);
    b_cov = cov(nData)+eps*eye(size(nData,2));
end
    
%Whiten Data
[U, D, V] = svd(b_cov);
sig_inv_half = D^(-1/2)*U';
dataBagsWhitened = {};
for i = 1:nBags
    m_minus = dataBags{i} - repmat(b_mu, [size(dataBags{i},1), 1]);
    m_scale = m_minus*sig_inv_half';
    if(parameters.methodFlag)
        denom = sqrt(repmat(sum(m_scale.*m_scale, 2), [1, nDim]));
        dataBagsWhitened{i} = m_scale./denom;
    else
        dataBagsWhitened{i} = m_scale;
    end
end

%Set up initial target signature, choose the best sample from the positive bags
pDataBags = dataBagsWhitened(labels ==  parameters.posLabel);
nDataBags = dataBagsWhitened(labels == parameters.negLabel);

if(parameters.initType == 1)
    [init_t, optObjVal, pBagsMax] = init1(pDataBags, nDataBags, parameters);
elseif(parameters.initType == 2)
    [init_t, optObjVal, pBagsMax] = init2(pDataBags, nDataBags, parameters);
else
    [init_t, optObjVal, pBagsMax] = init3(pDataBags, nDataBags, parameters);
end
optTarget = init_t;


%Precompute term 2 in update equation
nMean = zeros(length(nDataBags), nDim); 
for j = 1:length(nDataBags)  
    nData = nDataBags{j};
    nMean(j,:) = mean(nData)';
end
nMean = mean(nMean);

%Train target signature
objTracker(1).val = optObjVal;
objTracker(1).target = optTarget;
iter = 1;
continueFlag = 1; 
while(continueFlag && iter < parameters.maxIter  )
    iter = iter + 1;

    if(nPBags > 1)
        pMean = mean(pBagsMax);
    else
        pMean = pBagsMax;
    end
    t = pMean - nMean;
    optTarget = t/norm(t);

    %Update Objective and Determine the max points in each bag
    [optObjVal, pBagsMax] = evalObjectiveWhitened(pDataBags, nDataBags, optTarget,parameters.softmaxFlag);
    
    if(any([objTracker.val] == optObjVal))
        loc = ([objTracker.val] == optObjVal);
        if(length(loc) > 1)
            loc = loc(end);
        end
        if(~sum(abs(objTracker(loc).target - optTarget)))
            continueFlag = 0;
            disp(['Stopping at iter: ', num2str(iter)])
        end;
    end
    
    objTracker(iter).val = optObjVal;
    objTracker(iter).target = optTarget;

end

%Undo whitening
optTarget = (optTarget*D^(1/2)*V'); 
optTarget = optTarget/norm(optTarget);

init_t = (init_t*D^(1/2)*V'); 
init_t = init_t/norm(init_t);
% figure; plot(init_t);

end

function [objVal, pConfMax] = evalObjectiveWhitened(pDataBags, nDataBags, target, softMaxFlag)
numDim = size(pDataBags{1},2); 
pConfBags = zeros(length(pDataBags),2);
pConfMax = zeros(length(pDataBags), numDim);
for j = 1:length(pDataBags)  
    pData = pDataBags{j};
    pConf = sum(pData.*repmat(target, [size(pData,1), 1]),2);
    [maxConfs{1}(j), loc] = max(pConf); 
    if(softMaxFlag)
        w = exp(pConf)/sum(exp(pConf));
        pConfBags(j) = sum(pConf.*w); 
        pConfMax(j,:) = sum(repmat(w, [1, numDim]).*pData)';
    else
        pConfBags(j) = maxConfs{1}(j); 
        pConfMax(j,:) = pData(loc,:)';
    end
end

nConfBags = zeros(length(nDataBags),2);
for j = 1:length(nDataBags)  
    nData = nDataBags{j};
    nConf = sum(nData.*repmat(target, [size(nData,1), 1]),2);
    nConfBags(j) = mean(nConf);
end

objVal = mean(pConfBags(:,1)) - mean(nConfBags(:,1));

end

function [init_t, optObjVal, pBagsMax] = init1(pDataBags, nDataBags, parameters)

pData = vertcat(pDataBags{:});

tempObjVal = zeros(1,size(pData,1));
for j = 1:size(pData,1) %if large amount of data, can make this parfor loop
    optTarget = pData(j,:);
    tempObjVal(j) = evalObjectiveWhitened(pDataBags, nDataBags, optTarget,parameters.softmaxFlag);
end
[~, opt_loc] = max(tempObjVal);
optTarget = pData(opt_loc,:);
optTarget = optTarget/norm(optTarget);
init_t = optTarget;
[optObjVal, pBagsMax] = evalObjectiveWhitened(pDataBags, nDataBags, optTarget,parameters.softmaxFlag);

end


function [init_t, optObjVal, pBagsMax] = init2(pDataBags, nDataBags, parameters)

%Select data point with smallest cosine of the vector angle to background
%data
pData = vertcat(pDataBags{:});
nData = vertcat(nDataBags{:});
nDim = size(pData,2);
pdenom = sqrt(repmat(sum(pData.*pData, 2), [1, nDim]));
ndenom = sqrt(repmat(sum(nData.*nData, 2), [1, nDim]));
pData = pData./pdenom;
nData = nData./pdenom;

nDataMean = mean(nData);

tempObjVal = sum(pData.*repmat(nDataMean, [size(pData,1), 1]),2);
[~, opt_loc] = min(tempObjVal);
optTarget = pData(opt_loc,:);
optTarget = optTarget/norm(optTarget);
init_t = optTarget;
[optObjVal, pBagsMax] = evalObjectiveWhitened(pDataBags, nDataBags, optTarget,parameters.softmaxFlag);

end


function [init_t, optObjVal, pBagsMax] = init3(pDataBags, nDataBags, parameters)

pData = vertcat(pDataBags{:});
numP = size(pData,1);
rp = randperm(numP);
init_t = pData(rp(1),:);
[optObjVal, pBagsMax] = evalObjectiveWhitened(pDataBags, nDataBags, init_t,parameters.softmaxFlag);

end

