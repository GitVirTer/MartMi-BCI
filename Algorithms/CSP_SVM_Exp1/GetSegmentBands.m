function [RngMat, iCnt] = GetSegmentBands


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

wndWidth = 500;
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

wndWidth = 1000;
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

% wndWidth = 16;
% noverlap = wndWidth/2;
% for iHz = 4:noverlap:40
%     if(iHz+wndWidth < 40)
%         iCnt = iCnt + 1;
%         Rng = [iHz, iHz+wndWidth];
%         RngMat = [RngMat; Rng];
%     else
%         break;
%     end
% end
% 
% wndWidth = 32;
% noverlap = wndWidth/2;
% for iHz = 4:noverlap:40
%     if(iHz+wndWidth < 40)
%         iCnt = iCnt + 1;
%         Rng = [iHz, iHz+wndWidth];
%         RngMat = [RngMat; Rng];
%     else
%         break;
%     end
% end