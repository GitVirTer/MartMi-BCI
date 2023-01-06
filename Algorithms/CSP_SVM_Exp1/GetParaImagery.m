function ParaImagery = GetParaImagery(ImageData, ImageLabel, SelFlag, ParallelFlag, CSP_Config, nClass)

%% Filter
tic;
[ImageDataCell, ~, nBand] = FBNN_Filter(ImageData, []);
disp(['Filtering Completed, Time Comsumption:' num2str(toc) 's']);
%% Segment
tic;
[ImageDataSeg, ~, nSeg] = FBNN_Segment(ImageDataCell, []);
% clearvars DataTrainCell DataTestCell
disp(['Segmentation Completed, Time Comsumption:' num2str(toc) 's']);
%% FB特征提取
tic;
nCSP = 2;

if ParallelFlag
    parfor iSegBand = 1:nSeg*nBand
        iBand = mod(iSegBand,nBand);
        if iBand==0 iBand=nBand; end
        iSeg = floor((iSegBand-1)/nBand)+1;

        [VarMapImageData(iSegBand, :, :), ~, Wcsp{iSegBand}] = FilterBankFeatureExt(ImageDataSeg{iSeg}{iBand}, [], ImageLabel, [], nCSP, CSP_Config);
        disp(['Extrating Features... iSegBand = ' num2str(iSegBand)]);

    end
else
    for iSeg = 1:nSeg
        for iBand = 1:nBand
            iSegBand = (iSeg-1)*nBand+iBand;    %检查！！！
            [VarMapImageData(iSegBand, :, :), ~, Wcsp{iSegBand}] = FilterBankFeatureExt(ImageDataSeg{iSeg}{iBand}, [], ImageLabel, [], nCSP, CSP_Config);
            disp(['Extrating Features... iSegBand = ' num2str(iSegBand)]);
        end
    end    
end
% clearvars DataTrainSeg DataTestSeg

for iTrail = 1:size(VarMapImageData,3)
    VarMapData = VarMapImageData(:, :, iTrail);
%     VarMapData = VarMapData./max(max(VarMapData));
%     VarMapData = VarMapData(logical(PatCspIdx{iPat}), :);
    VarMapData = mapminmax(VarMapData,0,1);
%     FeatureMatTrain(:, :, 1, iTrail) = VarMapData;
    Feature2DMatTrain(:, :, 1, iTrail) = VarMapData;
    VarMapData = VarMapData';
    FeatureMatTrain(1, :, 1, iTrail) = VarMapData(:);
end

PatFeature.FeatureMatTrain = FeatureMatTrain;
PatFeature.Feature2DMatTrain = Feature2DMatTrain;
PatFeature.DataTrainLabel = ImageLabel;

ParaImagery.Wcsp = Wcsp;    % Wcsp

disp(['Feature Extraction Completed, Time Comsumption:' num2str(toc) 's']);

%% Load Feature Selection Index
accMat = {};
PatCspIdx = {};
if SelFlag
    tic;
    indices_train = crossvalind('Kfold',numel(find(double(ImageLabel)==1)),2);
    [DataSel.Train, DataSel.Test, DataSel.TrainLabel, DataSel.TestLabel] = GetFoldData(ImageData, ImageLabel, ~(indices_train == 1), (indices_train == 1), nClass);
    [~, ~, ~, ~, ~, PatFeature_Sel] = FBNN_Excute_CPIII(DataSel.Train,DataSel.TrainLabel,DataSel.Test,DataSel.TestLabel,ParallelFlag,CSP_Config,1);
    [accMat{1}, kappaMat{1}] = FBNN_Excute_FeatureSelection(PatFeature_Sel,ParallelFlag,1);
    PlotTFMap_All(accMat);  %画出图像
    PatCspIdx = GetDesignatedFeatureIndex(accMat, 150);
    disp(['Feature Extraction Completed, Time Comsumption:' num2str(toc) 's']);
end
ParaImagery.accMat = accMat;
ParaImagery.PatCspIdx = PatCspIdx;

for iTrail = 1:size(Feature2DMatTrain,4)
    if SelFlag    
        VarMapData = Feature2DMatTrain(logical(PatCspIdx{1}), :, 1, iTrail);
    else
        VarMapData = Feature2DMatTrain(:, :, 1, iTrail);
    end
    VarMapData = VarMapData';        
    FeatureMatTrainForSVM(iTrail,:) = VarMapData(:);
    
end

%% Train SVM
[~, ~, ~, ~, ~, Trained_LSVM] = trainSVM_2Class(FeatureMatTrainForSVM, [], ImageLabel, []);
ParaImagery.Trained_LSVM = Trained_LSVM;
