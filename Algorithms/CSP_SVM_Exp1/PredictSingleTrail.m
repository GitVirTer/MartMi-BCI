function PredictLabel = PredictSingleTrail(ParaImagery, ImageData, SelFlag, CSP_Config)
ParallelFlag = 0;
% CSP_Config.Wcsp = ParaImagery.Wcsp;
%% Filter
% tic;
[ImageDataCell, ~, nBand] = FBNN_Filter(ImageData, []);
% disp(['滤波完成 ，用时' num2str(toc) '秒']);
%% Segment
% tic;
[ImageDataSeg, ~, nSeg] = FBNN_Segment(ImageDataCell, []);
% clearvars DataTrainCell DataTestCell
% disp(['分段完成 ，用时' num2str(toc) '秒']);

nCSP = 2;

% if ParallelFlag
%     parfor iSegBand = 1:nSeg*nBand
%         iBand = mod(iSegBand,nBand);
%         if iBand==0 iBand=nBand; end
%         iSeg = floor((iSegBand-1)/nBand)+1;
%         CSP_Config.Wcsp = ParaImagery.Wcsp{iSegBand};
%         [VarMapImageData(iSegBand, :, :), ~, ~] = FilterBankFeatureExt(ImageDataSeg{iSeg}{iBand}, [], ImageLabel, [], nCSP, CSP_Config);
%         disp(['Extrating Features... iSegBand = ' num2str(iSegBand)]);
% 
%     end
% else
for iSeg = 1:nSeg
    for iBand = 1:nBand
        iSegBand = (iSeg-1)*nBand+iBand;    %检查！！！
        CSP_Config.Wcsp = ParaImagery.Wcsp{iSegBand};
        [VarMapImageData(iSegBand, :, :), ~, ~] = FilterBankFeatureExt(ImageDataSeg{iSeg}{iBand}, [], [], [], nCSP, CSP_Config);
%         disp(['Extrating Features... iSegBand = ' num2str(iSegBand)]);
    end
end    
% end

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

for iTrail = 1:size(Feature2DMatTrain,4)
    if SelFlag    
        VarMapData = Feature2DMatTrain(logical(ParaImagery.PatCspIdx{1}), :, 1, iTrail);
    else
        VarMapData = Feature2DMatTrain(:, :, 1, iTrail);
    end
    VarMapData = VarMapData';        
    FeatureMatTrainForSVM(iTrail,:) = VarMapData(:);
    
end

% [~, ~, ~, ~, ~, Trained_LSVM] = trainSVM_2Class(FeatureMatTrainForSVM, [], ImageLabel, []);
PredictLabel = predict(ParaImagery.Trained_LSVM, FeatureMatTrainForSVM);

