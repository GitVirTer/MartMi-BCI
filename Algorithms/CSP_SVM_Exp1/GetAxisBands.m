function [TimeAxis, FreqAxis, nTime, nFreq] = GetAxisBands
%% TimeBands
iCnt = 0;
RngMat = [];

wndWidth = 125;
noverlap = wndWidth;
for iHz = 1:noverlap:1000
    if(iHz+wndWidth-1 <= 1000)
        iCnt = iCnt + 1;
        Rng = [iHz, iHz+wndWidth-1];
        RngMat = [RngMat; Rng];
    else
        break;
    end
end

TimeAxis = RngMat;
nTime = iCnt;

%% FreqBands
iCnt = 0;
RngMat = [];

wndWidth = 2;
noverlap = wndWidth;
for iHz = 4:noverlap:36
    if(iHz+wndWidth-1 <= 1000)
        iCnt = iCnt + 1;
        Rng = [iHz, iHz+wndWidth-1];
        RngMat = [RngMat; Rng];
    else
        break;
    end
end

FreqAxis = RngMat;
nFreq = iCnt;
