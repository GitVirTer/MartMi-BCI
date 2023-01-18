function [net, Acc, Kappa, tp2, Wcsp, PatFeature] = FBNN_Excute_CPIII(DataTrain, DataTrainLabel, DataTest, DataTestLabel, ParallelFlag, CSP_Config, iPat)

net = {}; Acc = 0; Kappa = 0; tp2 = {};
% PatCspIdx = {};
% load PatCspIdx_CPIII_1.mat;
% load('PatCspIdx_AllPat_CPIII_First150.mat', 'PatCspIdx');

%% Filter
tic;
[DataTrainCell, DataTestCell, nBand] = FBNN_Filter(DataTrain, DataTest);
% save(['ButterFilteredData_Diff_Pat', num2str(iPat), '.mat'], 'DataTrainCell', 'DataTestCell', 'nBand', '-v7.3');

% clearvars DataTrain DataTest
% load(['ButterFilteredData_CPIII_Pat', num2str(iPat), '.mat'], 'DataTrainCell', 'DataTestCell', 'nBand');
% load(['ButterFilteredData_Pat', num2str(iPat), '.mat'], 'DataTrainCell', 'DataTestCell', 'nBand');
% load(['ButterFilteredData_CPIII_Pat', num2str(iPat), '.mat'], 'DataTrainCell', 'DataTestCell', 'nBand');
toc;
%% Segment
tic;
% for i = 1:numel(DataTrainCell)
%     DataTrainCell{i} = DataTrainCell{i}(:,1:22,:);
% end
% for i = 1:numel(DataTestCell)
%     DataTestCell{i} = DataTestCell{i}(:,1:22,:);
% end
[DataTrainSeg, DataTestSeg, nSeg] = FBNN_Segment(DataTrainCell, DataTestCell);
% clearvars DataTrainCell DataTestCell
toc;
%% FB特征提取
tic;
nCSP = 2;

if ParallelFlag
    parfor iSegBand = 1:nSeg*nBand
        iBand = mod(iSegBand,nBand);
        if iBand==0 iBand=nBand; end
        iSeg = floor((iSegBand-1)/nBand)+1;

        [VarMapTrain(iSegBand, :, :), VarMapTest(iSegBand, :, :), Wcsp{iSegBand}] = FilterBankFeatureExt(DataTrainSeg{iSeg}{iBand}, DataTestSeg{iSeg}{iBand}, DataTrainLabel, DataTestLabel, nCSP, CSP_Config);
        disp(['Extracting Features... iSegBand = ' num2str(iSegBand)]);

    end
else
    for iSegBand = 1:nSeg*nBand
        iBand = mod(iSegBand,nBand);
        if iBand==0 iBand=nBand; end
        iSeg = floor((iSegBand-1)/nBand)+1;

        [VarMapTrain(iSegBand, :, :), VarMapTest(iSegBand, :, :), Wcsp{iSegBand}] = FilterBankFeatureExt(DataTrainSeg{iSeg}{iBand}, DataTestSeg{iSeg}{iBand}, DataTrainLabel, DataTestLabel, nCSP, CSP_Config);
        disp(['Extracting Features... iSegBand = ' num2str(iSegBand)]);

    end 
end
% clearvars DataTrainSeg DataTestSeg

for iTrail = 1:size(VarMapTrain,3)
    VarMapData = VarMapTrain(:, :, iTrail);
%     VarMapData = VarMapData./max(max(VarMapData));
%     VarMapData = VarMapData(logical(PatCspIdx{iPat}), :);
    VarMapData = mapminmax(VarMapData,0,1);
%     FeatureMatTrain(:, :, 1, iTrail) = VarMapData;
    Feature2DMatTrain(:, :, 1, iTrail) = VarMapData;
    VarMapData = VarMapData';
    FeatureMatTrain(1, :, 1, iTrail) = VarMapData(:);
end

for iTrail = 1:size(VarMapTest,3)
    VarMapData = VarMapTest(:, :, iTrail);
%     VarMapData = VarMapData./max(max(VarMapData));
%     VarMapData = VarMapData(logical(PatCspIdx{iPat}), :);
    VarMapData = mapminmax(VarMapData,0,1);
%     FeatureMatTest(:, :, 1, iTrail) = VarMapData;
    Feature2DMatTest(:, :, 1, iTrail) = VarMapData;
    VarMapData = VarMapData';
    FeatureMatTest(1, :, 1, iTrail) = VarMapData(:);
    
end

PatFeature.FeatureMatTrain = FeatureMatTrain;
PatFeature.Feature2DMatTrain = Feature2DMatTrain;
PatFeature.DataTrainLabel = DataTrainLabel;

PatFeature.FeatureMatTest = FeatureMatTest;
PatFeature.Feature2DMatTest = Feature2DMatTest;
PatFeature.DataTestLabel = DataTestLabel;

% DataTrain_Format = FeatureMatTrain;
% DataTest_Format = FeatureMatTest;
toc;
%% 训练

% [net, sumh2, suma2] = trainNet2_ButterFB(DataTrain_Format, DataTest_Format, DataTrainLabel, DataTestLabel);

%% Confusion Matrix and Kappa Coefficient

% tp2 = zeros(4,4);
% for numx = 1:length(sumh2)
%     tp2(suma2(numx),sumh2(numx)) = tp2(suma2(numx),sumh2(numx)) + 1;
% end
% 
% a = sum(tp2);
% b = sum(tp2.');
% N = sum(sum(tp2));
% accSum = sum(diag(tp2));
% 
% p0 = accSum/N;
% pe = a*b'/(N^2);
% 
% Kappa = (p0-pe)/(1-pe);
% Acc = sum(diag(tp2))/sum(sum(tp2));
% 
% disp(['Kappa Coefficient is: ' num2str(Kappa)]);
