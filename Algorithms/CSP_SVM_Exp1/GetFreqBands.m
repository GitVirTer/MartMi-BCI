function [RngMat, iCnt] = GetFreqBands


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

wndWidth = 8;
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

wndWidth = 16;
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

wndWidth = 32;
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