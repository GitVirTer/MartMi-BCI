function [DataTrainSeg, DataTestSeg, nSeg] = WavPairCSPSeg_Segment(DataTrain, DataTest, SegLen, noverlap)
iSeg = 0;
for iDot = 1:noverlap:size(DataTrain{1},3)-SegLen+1
    iSeg = iSeg + 1;
    for iBand = 1:numel(DataTrain)
        DataTrainSeg{iSeg}{iBand} = DataTrain{iBand}(:, :, iDot:iDot+SegLen-1);
    end
end

if ~isempty(DataTest)
    iSeg = 0;
    for iDot = 1:noverlap:size(DataTest{1},3)-SegLen+1
        iSeg = iSeg + 1;
        for iBand = 1:numel(DataTest)
            DataTestSeg{iSeg}{iBand} = DataTest{iBand}(:, :, iDot:iDot+SegLen-1);
        end
    end
else
    DataTestSeg = [];
end

if iDot+SegLen-1 < size(DataTrain{1},3)
    warning(['Missed ' num2str(size(DataTrain{1},3)-(iDot+SegLen-1)) ' Dot!']);
end

nSeg = iSeg;