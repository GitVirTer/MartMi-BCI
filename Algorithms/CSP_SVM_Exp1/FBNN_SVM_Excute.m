function [Acc, Kappa, tp2] = FBNN_SVM_Excute(PatFeature,accMat,SelFlag,iPat)
%% Load Feature Selection Index
% PatCspIdx = {};
% load PatCspIdx_CPIII_1.mat;
% load('PatCspIdx_AllPat_CPIII_First150.mat', 'PatCspIdx');
if SelFlag
    PatCspIdx = GetDesignatedFeatureIndex(accMat, 150);
end
% PatCspIdx = GetFeatureIndex;
% load('PatFeature_22Ch.mat', 'PatFeature');
% load('PatFeature_divCSP_CPIII_IVa.mat', 'PatFeature');
% FeatureMatTrain = PatFeature.FeatureMatTrain;
Feature2DMatTrain = PatFeature.Feature2DMatTrain;
DataTrainLabel = PatFeature.DataTrainLabel;

% FeatureMatTest = PatFeature.FeatureMatTest;
Feature2DMatTest = PatFeature.Feature2DMatTest;
DataTestLabel = PatFeature.DataTestLabel;

%% Filter
% % [DataTrainCell, DataTestCell, nBand] = FBNN_Filter(DataTrain, DataTest);
% % save(['ButterFilteredData_Diff_Pat', num2str(iPat), '.mat'], 'DataTrainCell', 'DataTestCell', 'nBand', '-v7.3');
% 
% clearvars DataTrain DataTest
% % load(['ButterFilteredData_CPIII_Pat', num2str(iPat), '.mat'], 'DataTrainCell', 'DataTestCell', 'nBand');
% % load(['ButterFilteredData_Pat', num2str(iPat), '.mat'], 'DataTrainCell', 'DataTestCell', 'nBand');
% load(['ButterFilteredData_Pat', num2str(iPat), '.mat'], 'DataTrainCell', 'DataTestCell', 'nBand');

%% Segment
% for i = 1:numel(DataTrainCell)
%     DataTrainCell{i} = DataTrainCell{i}(:,1:22,:);
% end
% for i = 1:numel(DataTestCell)
%     DataTestCell{i} = DataTestCell{i}(:,1:22,:);
% end
% [DataTrainSeg, DataTestSeg, nSeg] = FBNN_Segment(DataTrainCell, DataTestCell);
% clearvars DataTrainCell DataTestCell

%% FB特征提取
tic;
% nCSP = 2;
% for iSeg = 1:nSeg
%     for iBand = 1:nBand
%         iSegBand = (iSeg-1)*nBand+iBand;    %检查！！！
%         [VarMapTrain(iSegBand, :, :), VarMapTest(iSegBand, :, :), Wcsp{iSegBand}] = FilterBankFeatureExt(DataTrainSeg{iSeg}{iBand}, DataTestSeg{iSeg}{iBand}, DataTrainLabel, DataTestLabel, nCSP);
%     end
% end
% clearvars DataTrainSeg DataTestSeg
% 
% for iTrail = 1:size(VarMapTrain,3)
%     VarMapData = VarMapTrain(:, :, iTrail);
% %     VarMapData = VarMapData./max(max(VarMapData));
% %     VarMapData = VarMapData(logical(PatCspIdx{iPat}), :);
%     VarMapData = mapminmax(VarMapData,0,1);
%     Feature2DMatTrain(:, :, 1, iTrail) = VarMapData;
%     VarMapData = VarMapData';
%     
% %     FeatureMatTrain(:, :, 1, iTrail) = VarMapData;
%     FeatureMatTrain(1, :, 1, iTrail) = VarMapData(:);
%     
%     VarMapData = VarMapTest(:, :, iTrail);
% %     VarMapData = VarMapData./max(max(VarMapData));
% %     VarMapData = VarMapData(logical(PatCspIdx{iPat}), :);
%     VarMapData = mapminmax(VarMapData,0,1);
%     Feature2DMatTest(:, :, 1, iTrail) = VarMapData;
%     VarMapData = VarMapData';
%     
% %     FeatureMatTest(:, :, 1, iTrail) = VarMapData;
%     FeatureMatTest(1, :, 1, iTrail) = VarMapData(:);
%     
% end

% FeatureMatTrain = PatFeature{iPat}.FeatureMatTrain;
% Feature2DMatTrain = PatFeature{iPat}.Feature2DMatTrain;
for iTrail = 1:size(Feature2DMatTrain,4)
    if SelFlag    
        VarMapData = Feature2DMatTrain(logical(PatCspIdx{iPat}), :, 1, iTrail);
    else
        VarMapData = Feature2DMatTrain(:, :, 1, iTrail);
    end
%     VarMapData = diff(VarMapData,[],2);
    VarMapData = VarMapData';
%     sumMat = [];
%     for i = 1:6
%         subMat = VarMapData((i-1)*4+1:i*4,:);
%         subMat = mapminmax(subMat,0,1);
%         sumMat = cat(1,sumMat,subMat);
%     end
%     VarMapData = sumMat;
%     VarMapData = mapminmax(VarMapData,0,1);

%     FeatureMatTrain(1, :, 1, iTrail) = mapminmax(VarMapData(:),0,1);
%     FeatureMatTrain(iTrail,:) = mapminmax((diff(VarMapData(:))),0,1);
    
        
    FeatureMatTrain(iTrail,:) = VarMapData(:);
    
end

% % Attractor Metagene Feature Selection
% A = reshape(FeatureMatTrain,size(FeatureMatTrain,2),size(FeatureMatTrain,4));
% [~,W] = metafeatures(A);
% [~, SortIdx] = sort(W, 'descend');
% FeatureMatTrain = FeatureMatTrain(:,SortIdx(1:600),:,:);

% FeatureMatTest = PatFeature{iPat}.FeatureMatTest;
% Feature2DMatTest = PatFeature{iPat}.Feature2DMatTest;
for iTrail = 1:size(Feature2DMatTest,4)
    if SelFlag    
        VarMapData = Feature2DMatTest(logical(PatCspIdx{iPat}), :, 1, iTrail);
    else
        VarMapData = Feature2DMatTest(:, :, 1, iTrail);
    end
%     VarMapData = diff(VarMapData,[],2);
    VarMapData = VarMapData';
%     sumMat = [];
%     for i = 1:6
%         subMat = VarMapData((i-1)*4+1:i*4,:);
%         subMat = mapminmax(subMat,0,1);
%         sumMat = cat(1,sumMat,subMat);
%     end
%     VarMapData = sumMat;
%     VarMapData = mapminmax(VarMapData,0,1);
    
%     FeatureMatTest(1, :, 1, iTrail) = mapminmax(VarMapData(:),0,1);
%     FeatureMatTest(iTrail,:) = mapminmax((diff(VarMapData(:))),0,1);
    FeatureMatTest(iTrail,:) = VarMapData(:);
    
end

% % Attractor Metagene Feature Selection
% FeatureMatTest = FeatureMatTest(:,SortIdx(1:600),:,:);

DataTrain_Format = FeatureMatTrain;
DataTest_Format = FeatureMatTest;
toc;

%% 训练
% [B,FitInfo] = lasso(DataTrain_Format, DataTrainLabel, 'Lambda',0.0101, 'Alpha', 0.5);
% for i = 1:size(B,2) data = B(:,i); Bnot0(i) = numel(data(data~=0)); end
% disp(['Not zero number: ' num2str(Bnot0)]);

[~, ~, ~, PredScore4, numB] = trainSVM_2Class(DataTrain_Format, DataTest_Format, DataTrainLabel, DataTestLabel);

for iB = 1:numB
%     PredScore = [PredScore1{iB}(:,2) PredScore2{iB}(:,2) PredScore3{iB}(:,2) PredScore4{iB}(:,2)];
    PredScore = PredScore4{iB};
    [~, sumh2] = max(PredScore,[],2);
    suma2 = DataTestLabel;

%% Confusion Matrix and Kappa Coefficient

    tp2{iB} = zeros(2,2);
    for numx = 1:length(sumh2)
        tp2{iB}(suma2(numx),sumh2(numx)) = tp2{iB}(suma2(numx),sumh2(numx)) + 1;
    end

    [Kappa(iB), Acc(iB)] = GetKappaAcc(tp2{iB});

    disp(['Acc is: ' num2str(Acc(iB))]);
end
