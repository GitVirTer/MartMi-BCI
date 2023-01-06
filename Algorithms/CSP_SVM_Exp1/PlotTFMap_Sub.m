function PlotTFMap_Sub(kappaMat)

[AxisTime, AxisFreq, nTime, nFreq] = GetAxisBands;
[RngMatSeg, nBandSeg] = GetSegmentBandsTest;
[RngMatFreq, nBandFreq] = GetFreqBandsTest;

for iTime = 1:nTime
    Time = AxisTime(iTime, :);
    thr = (Time(1)>=RngMatSeg(:,1))&(Time(2)<=RngMatSeg(:,2));
    IdxTime = find(thr == 1);
    TimeCell{iTime} = IdxTime;
end

for iFreq = 1:nFreq
    Freq = AxisFreq(iFreq, :);
    thr = (Freq(1)>=RngMatFreq(:,1))&(Freq(2)<=RngMatFreq(:,2));
    IdxFreq = find(thr == 1);
    FreqCell{iFreq} = IdxFreq;
end

%% Start Fuse Matrices
% load('kappaMap_60.mat', 'kappaMat')
% kappaMat = accMat;
% TMat = zeros(numel(FreqCell), numel(TimeCell));
% FMat = zeros(numel(FreqCell), numel(TimeCell));
TFMat = zeros(numel(FreqCell), numel(TimeCell));
for iPat = 1:numel(kappaMat)
    kappaMatPat = kappaMat{iPat}';
    kappaMatPat = kappaMatPat(1:16, 1:7);
    for iT = 1:numel(TimeCell)
        TMat{iPat}(:, iT) = sum((1/numel(TimeCell{iT})).*kappaMatPat(:, TimeCell{iT}),2);
    end
end

for iPat = 1:numel(kappaMat)
    TMatPat = TMat{iPat};
    for iF = 1:numel(FreqCell)
        FMat{iPat}(iF, :) = sum((1/numel(FreqCell{iF})).*TMatPat(FreqCell{iF}, :),1);
    end
end

TFMat = FMat;

figure;
for i = 1:numel(TFMat)
    subplot(1,numel(TFMat),i)
    imagesc(0:0.5:4, 4:2:36, TFMat{i});set(gca,'ydir','normal');
end


function [RngMat, iCnt] = GetSegmentBandsTest


iCnt = 0;
RngMat = [];

wndWidth = 250;
noverlap = wndWidth/2;
for iHz = 1:noverlap:1000
    if(iHz+wndWidth-1 <= 1000)
        iCnt = iCnt + 1;
        Rng = [iHz, iHz+wndWidth-1];
        RngMat = [RngMat; Rng];
    else
        break;
    end
end

end

function [RngMat, iCnt] = GetFreqBandsTest

iCnt = 0;
RngMat = [];

wndWidth = 4;
noverlap = wndWidth/2;
for iHz = 4:noverlap:40
    if(iHz+wndWidth < 40)
        iCnt = iCnt + 1;
        Rng = [iHz, iHz+wndWidth];
        RngMat = [RngMat; Rng];
    else
        break;
    end
end
end

end