function PatCspIdx = GetDesignatedFeatureIndex(kappaMat, NumFeature)
% load('FBNN_Results_nCSP2_2.mat');
% load('FeatureSelection_TrainResMat_AllPat.mat', 'kappaMat')
% load('kappaMap_divCSP_36.mat', 'kappaMat')
% load('divCSP_CPIV_Train_Half_KappaMat.mat', 'kappaMat')
% load('divCSP_CPIII_Half_KappaMat.mat', 'kappaMat')

for iPat = 1:numel(kappaMat)
    
    
    nSeg = size(kappaMat{iPat},1);
    nBand = size(kappaMat{iPat},2);
    for iSeg = 1:nSeg
        for iBand = 1:nBand
            iFeature = (iSeg-1)*nBand+iBand;
            kappaVector{iPat}(iFeature) = kappaMat{iPat}(iSeg, iBand);
        end
    end
    
%     kappaVector{iPat}(iFeature) = smooth(kappaVector{iPat}(iFeature), 5);       % ¼ÓÆ½»¬ÂË²¨£¡£¡£¡
    
    segIdx = zeros(1, nSeg*nBand);
    
    [~, Idx] = sort(kappaVector{iPat}, 'descend');
    segIdx(Idx(1:NumFeature)) = 1;
    
    PatCspIdx{iPat} = segIdx;
    
%     subplot(9,1,iPat)
% %     WMapPlot(:,:,1) = WMap';
% %     WMapPlot(:,:,2) = WMap';
% %     WMapPlot(:,:,3) = WMap';
% %     image(WMapPlot);
%     plot(WMapVector);
end