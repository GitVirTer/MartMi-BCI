function [accMat, kappaMat] = FBNN_Excute_FeatureSelection(PatFeature, ParallelFlag, iPat)
%% Transfer Data
% FeatureMatTrain = PatFeature.FeatureMatTrain;
Feature2DMatTrain = PatFeature.Feature2DMatTrain;
DataTrainLabel = PatFeature.DataTrainLabel;

% FeatureMatTest = PatFeature.FeatureMatTest;
Feature2DMatTest = PatFeature.Feature2DMatTest;
DataTestLabel = PatFeature.DataTestLabel;

VarMapTrain = reshape(Feature2DMatTrain, size(Feature2DMatTrain,1), size(Feature2DMatTrain,2), size(Feature2DMatTrain,4));
VarMapTest = reshape(Feature2DMatTest, size(Feature2DMatTest,1), size(Feature2DMatTest,2), size(Feature2DMatTest,4));
%% Normalization Data
for iTrail = 1:size(VarMapTrain,3)
    VarMapData = VarMapTrain(:, :, iTrail);
%     VarMapData = VarMapData';
%     VarMapData = mapminmax(VarMapData,0,1);
%     VarMapData = VarMapData';
    FeatureMatTrain(:, :, 1, iTrail) = VarMapData;

end

for iTrail = 1:size(VarMapTest,3)
    VarMapData = VarMapTest(:, :, iTrail);
%     VarMapData = VarMapData';
%     VarMapData = mapminmax(VarMapData,0,1);
%     VarMapData = VarMapData';
    FeatureMatTest(:, :, 1, iTrail) = VarMapData;
end

DataTrain_Format = FeatureMatTrain;
DataTest_Format = FeatureMatTest;
toc;

%% Feature Seletion
% ParallelFlag = 1;
nSeg = 11; 
nBand = 27;
accMat = zeros(nSeg, nBand);
kappaMat = zeros(nSeg, nBand);
if ParallelFlag
%     [accMat, kappaMat] = Parallel_trainSingleFeature(DataTrain_Format,DataTrainLabel,DataTest_Format,DataTestLabel);
    parfor iFeature = 1:nSeg*nBand
        iBand = mod(iFeature,nBand);
        if iBand==0 iBand=nBand; end
        iSeg = floor((iFeature-1)/nBand)+1;

        FeatureTrain = DataTrain_Format(iFeature, :, :, :);
        FeatureTest = DataTest_Format(iFeature, :, :, :);
        %% ÑµÁ·
        [net_par, sumh2, suma2] = trainNet2_ButterFeature(FeatureTrain, FeatureTest, DataTrainLabel, DataTestLabel);

        %% Confusion Matrix and Kappa Coefficient
        tp2_par = zeros(4,4);
        for numx = 1:length(sumh2)
            tp2_par(suma2(numx),sumh2(numx)) = tp2_par(suma2(numx),sumh2(numx)) + 1;
        end

        [Kappa_par, Acc_par] = GetKappaAcc(tp2_par);

        accMatVec(iFeature) = Acc_par;
        kappaMatVec(iFeature) = Kappa_par;    
    %     accMat(iSeg, iBand) = Acc;
    %     kappaMat(iSeg, iBand) = Kappa;

        disp(['Pat: ' num2str(iPat) ': ' num2str(iFeature) '/' num2str(nSeg*nBand) ' , Accuracy is: ' num2str(Acc_par) ' , Kappa Coefficient is: ' num2str(Kappa_par)]);

    %     end
    end

    accMat = zeros(nSeg, nBand);
    kappaMat = zeros(nSeg, nBand);
    for iFeature = 1:nSeg*nBand
        iBand = mod(iFeature,nBand);
        if iBand==0 iBand=nBand; end
        iSeg = floor((iFeature-1)/nBand)+1;

        accMat(iSeg, iBand) = accMatVec(iFeature);
        kappaMat(iSeg, iBand) = kappaMatVec(iFeature);
    end    

else
    for iSeg = 1:nSeg
        for iBand = 1:nBand
        iFeature = (iSeg-1)*nBand+iBand;

        FeatureTrain = DataTrain_Format(iFeature, :, :, :);
        FeatureTest = DataTest_Format(iFeature, :, :, :);
        %% ÑµÁ·
        [net, sumh2, suma2] = trainNet2_ButterFeature(FeatureTrain, FeatureTest, DataTrainLabel, DataTestLabel);

        %% Confusion Matrix and Kappa Coefficient
        tp2 = zeros(4,4);
        for numx = 1:length(sumh2)
            tp2(suma2(numx),sumh2(numx)) = tp2(suma2(numx),sumh2(numx)) + 1;
        end

        [Kappa, Acc] = GetKappaAcc(tp2);

%         accMatVec(iFeature) = Acc;
%         kappaMatVec(iFeature) = Kappa;    
        accMat(iSeg, iBand) = Acc;
        kappaMat(iSeg, iBand) = Kappa;

        disp(['Pat: ' num2str(iPat) ': ' num2str(iFeature) '/' num2str(nSeg*nBand) ' , Accuracy is: ' num2str(Acc) ' , Kappa Coefficient is: ' num2str(Kappa)]);

        end
    end    
    
end


