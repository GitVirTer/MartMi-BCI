function [DataTrainCell, DataTestCell, nBand] = FBNN_Filter(DataTrain, DataTest)
% tic;
% DataTrain, DataTest: 1x* cell
N = 5;
reSAMP = 250;
[RngMat, nBand] = GetFreqBands;
for iBand = 1:nBand
    [DataTrainCell{iBand}, DataTestCell{iBand}] = FilterBankProc(DataTrain, DataTest, RngMat(iBand, :), N, reSAMP);
end
% [DataTrainCell, DataTestCell] = WavPairCSPSeg_Filter(DataTrain, DataTest);
% nBand = numel(DataTrainCell);

% toc;
% disp('Filtering Completed!');