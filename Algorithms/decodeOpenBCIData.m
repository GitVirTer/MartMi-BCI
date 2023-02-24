function [UnpackedData, UnpackedDataRaw] = decodeOpenBCIData(ObjSerial, BufferSize)
%This function is designed for decoding OpenBCI Data.
%If you want to use other EEG recording device, please modify this function
%to obtain 'UnpackedData' and 'UnpackedDataRaw'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PackageLen = 33;
nCh = 8;
scale_fac_uVolts_per_count = 0.022351744455307063;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read Data From Buffer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recdta = fread(ObjSerial, BufferSize, 'uchar');   % Receive serial data from OpenBCI (33 Bytes/Package, 256 Package/Second)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unpack Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A0 = find(recdta==hex2dec('A0'));
iPackage = 0;
for iSt = 1:numel(A0)
    if (A0(iSt)+PackageLen-1<=numel(recdta))&&(recdta(A0(iSt)+PackageLen-1)==192)    % 'C0'=192
        iPackage = iPackage+1;
        PackArr(iPackage,:) = recdta(A0(iSt):A0(iSt)+PackageLen-1);
    end
end
PackArr = PackArr(:,3:26);
UnpackedData = zeros(nCh,size(PackArr,1));
for iCh = 1:nCh
    UnpackedData(iCh,:) = UnpackData(PackArr(:,(iCh-1)*3+1:iCh*3));
end
UnpackedData = UnpackedData.*scale_fac_uVolts_per_count;
UnpackedDataRaw = UnpackedData;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
