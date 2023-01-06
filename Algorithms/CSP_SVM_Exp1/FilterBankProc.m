function [DataTrainTrans, DataTestTrans] = FilterBankProc(DataTrain, DataTest, FreqRng, N, reSAMP)
DataTrainTrans = zeros(size(DataTrain));
for iTrail = 1:size(DataTrain, 1)
    for iCh = 1:size(DataTrain, 2)
        DataTrainTrans(iTrail, iCh, :) = ButterFilter_Mat(reshape(DataTrain(iTrail, iCh, :), 1, numel(DataTrain(iTrail, iCh, :))), FreqRng, N, reSAMP);
    end
end

if ~isempty(DataTest)
    DataTestTrans = zeros(size(DataTest));
    for iTrail = 1:size(DataTest, 1)
        for iCh = 1:size(DataTest, 2)
            DataTestTrans(iTrail, iCh, :) = ButterFilter_Mat(reshape(DataTest(iTrail, iCh, :), 1, numel(DataTest(iTrail, iCh, :))), FreqRng, N, reSAMP);
        end
    end
    
else
    DataTestTrans = [];
end