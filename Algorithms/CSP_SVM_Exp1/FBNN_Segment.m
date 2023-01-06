function [DataTrainSeg, DataTestSeg, nSeg] = FBNN_Segment(DataTrainCell, DataTestCell)
% tic;
SegLen = 250;
noverlap = 125;
[DataTrainSeg1, DataTestSeg1, nSeg1] = WavPairCSPSeg_Segment(DataTrainCell, DataTestCell, SegLen, noverlap);

SegLen = 500;
noverlap = 250;
[DataTrainSeg2, DataTestSeg2, nSeg2] = WavPairCSPSeg_Segment(DataTrainCell, DataTestCell, SegLen, noverlap);

SegLen = 1000;
noverlap = 1000;
[DataTrainSeg3, DataTestSeg3, nSeg3] = WavPairCSPSeg_Segment(DataTrainCell, DataTestCell, SegLen, noverlap);

DataTrainSeg = cat(2, DataTrainSeg1, DataTrainSeg2, DataTrainSeg3);

if ~isempty(DataTestCell)
    DataTestSeg = cat(2, DataTestSeg1, DataTestSeg2, DataTestSeg3);
else
    DataTestSeg = [];
end

nSeg = nSeg1 + nSeg2 + nSeg3;

% toc;
% disp('Segmentation Completed!');