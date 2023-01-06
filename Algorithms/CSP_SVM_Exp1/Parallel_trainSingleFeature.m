function [accMat, kappaMat] = Parallel_trainSingleFeature(DataTrain_Format,DataTrainLabel,DataTest_Format,DataTestLabel)
nSeg = 11; 
nBand = 27;
accMat = zeros(nSeg, nBand);
kappaMat = zeros(nSeg, nBand);

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

    disp(['Pat: ' num2str(iPat) ': ' num2str(iFeature) '/' num2str(nSeg*nBand) ' , Accuracy is: ' num2str(Acc) ' , Kappa Coefficient is: ' num2str(Kappa)]);

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
