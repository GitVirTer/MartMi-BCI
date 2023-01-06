function varargout = RTEEGAP(varargin)
% RTEEGAP MATLAB code for RTEEGAP.fig
%      RTEEGAP, by itself, creates a new RTEEGAP or raises the existing
%      singleton*.
%
%      H = RTEEGAP returns the handle to a new RTEEGAP or the handle to
%      the existing singleton*.
%
%      RTEEGAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RTEEGAP.M with the given input arguments.
%
%      RTEEGAP('Property','Value',...) creates a new RTEEGAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RTEEGAP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RTEEGAP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RTEEGAP

% Last Modified by GUIDE v2.5 09-May-2022 13:00:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RTEEGAP_OpeningFcn, ...
                   'gui_OutputFcn',  @RTEEGAP_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before RTEEGAP is made visible.
function RTEEGAP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RTEEGAP (see VARARGIN)

global SerialPortNum;
global BaudRate;
global Filter;
global DrawMap;
global ZoomCoeff;
global SampleRate;
global DispSec;
global SecCnt;
global BufferSize
% global CircMaskData
global ElectrodeMap
global FFTMap
global PowerMap
global IntervalTime
global NumTryStartCom
global RecordData
global MotorImagery

CurrentPath = cd;
ThisToolRootPath = [CurrentPath, '\'];
addpath(genpath(ThisToolRootPath));
set(handles.EditImageFilePath, 'String', [CurrentPath, '\Data\MotorImageryData\']);
set(handles.EditImageConfigureFilePath, 'String', [CurrentPath, '\Data\Configure\Configure.mat']);

MotorImagery.RecordFlag = 0;
MotorImagery.RecordData = [];
MotorImagery.RecordRawData = [];
MotorImagery.EpochCnt = 0;
MotorImagery.PrepareTime = str2double(get(handles.EditPrepareTime, 'String'));              % 准备时间
MotorImagery.ImageLastTime = str2double(get(handles.EditImageLastTime, 'String'));          % 持续时间
MotorImagery.CueTime = str2double(get(handles.EditCueTime, 'String'));                      % 提示时间
MotorImagery.AllTime = str2double(get(handles.EditAllTime, 'String'));                      % 末尾时间
MotorImagery.StageFlag = 0;                                                                % 时间段标志位
MotorImagery.Sample.Cnt = 0;                                                                    % 样本点数
MotorImagery.Sample.Stage1 = 0;
MotorImagery.Sample.Stage2 = 0;
MotorImagery.Sample.Stage3 = 0;
MotorImagery.Sample.Stage4 = 0;
MotorImagery.ImageFilePath = get(handles.EditImageFilePath, 'String');                      % 运动想象数据保存路径
MotorImagery.CurrentPath = CurrentPath;
MotorImagery.ImageSound.FilePath = [MotorImagery.CurrentPath, '\Data\Sound\cue.mp3'];
[MotorImagery.ImageSound.SoundCue, MotorImagery.ImageSound.SoundFs] = audioread(MotorImagery.ImageSound.FilePath);
MotorImagery.ImageConfigureFilePath = get(handles.EditImageConfigureFilePath, 'String');    % 配置文件路径
MotorImagery.ClassName = get(handles.ListBoxImageEvent, 'String');
MotorImagery.numClass = numel(MotorImagery.ClassName);
MotorImagery.ClassEpoch = str2double(get(handles.EditClassEpoch, 'String'));
MotorImagery.AllEpoch = MotorImagery.numClass*MotorImagery.ClassEpoch;                % 总试验次�?
MotorImagery.RemainEpoch = MotorImagery.AllEpoch;                                           % 剩余次数
MotorImagery.State = 'NOT Start';
MotorImagery.ParaImagery = [];                                              % 预测参数
set(handles.TextImageExperimentState, 'String', MotorImagery.State);    % 用此句实时更新状�?
ClassOrder = repmat(1:MotorImagery.numClass,[1 MotorImagery.ClassEpoch]);
ClassOrder = ClassOrder(randperm(numel(ClassOrder)));
MotorImagery.ClassOrder = ClassOrder;
MotorImagery.CurrentCnt = 0;
% MotorImagery.CurrentClass = MotorImagery.ClassOrder(1);
for iClass = 1:MotorImagery.numClass
    MotorImagery.ClassImage{iClass} = imread([CurrentPath, '\Data\Picture\' MotorImagery.ClassName{iClass} '.jpg']);
end
MotorImagery.PrepareImage = imread([CurrentPath, '\Data\Picture\Prepare.jpg']);
MotorImagery.NoCueImage = imread([CurrentPath, '\Data\Picture\NoCue.jpg']);
MotorImagery.UseBigCue = 1;
sound(MotorImagery.ImageSound.SoundCue, MotorImagery.ImageSound.SoundFs);   % 播放音乐
set(handles.CheckBoxUseBigCue, 'Value', MotorImagery.UseBigCue);
save(get(handles.EditImageConfigureFilePath, 'String'), 'MotorImagery');

image(handles.AxesImageSmallCue, MotorImagery.PrepareImage);
axis(handles.AxesImageSmallCue, 'off');
% image(handles.AxesImageBigCue, MotorImagery.PrepareImage);
plot(handles.AxesImageBigCue,0,0);
set(handles.AxesImageBigCue, 'Visible', 'off');
set(handles.uipanel1, 'Visible', 'on');
set(handles.uipanel2, 'Visible', 'on');
set(handles.uipanel3, 'Visible', 'on');

RecordData.NowData = [];
RecordData.RecFlag = 0;
RecordData.EventFlag = 0;
RecordData.Path = get(handles.EditSavePath,'String');
RecordData.Time = 0;
RecordData.Sample = 0;
RecordData.EventTime = [];
RecordData.EventSample = [];
RecordData.EventDescribe = {};

NumTryStartCom = 0;
SecCnt = 0;

% SerialPortNum = 'COM4';
SerialPortNum = get(handles.EditSerialPort,'String');
BaudRate = 115200;

SampleRate = 250;   % The sample rate of EEG stream
DispSec = 10;   % The length of EEG stream displayed
IntervalTime.DispMs = 500;   % The number of millisecond in refresh interval of displaying
IntervalTime.FFTMs = 1000;

nCh = 8;    % The number of channel
Fs  = SampleRate;  % Sampling Frequency
N   = 8;  % Order
Fc1 = 1;   % First Cutoff Frequency
Fc2 = 35;  % Second Cutoff Frequency

% IntervalTime.PSDMs = 1000;

BufferSize = round(33*SampleRate*(IntervalTime.DispMs/1000));    % 设置缓存区大�?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 设置脑电地形�? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BG_color = 240;
r=64;
m=2*r; %矩阵的行�?
n=2*r; %矩阵的列�?
% r=16;   %生成圆的半径
m1=-m/2:m/2-1;   %把圆心变到矩阵的中间
n1=-n/2:n/2-1;
[x,y]=meshgrid(m1,n1);
circle=x.^2+y.^2;   %计算出每�?点到圆心的距离的平方

% Gaussian Kernel
sig = 25;
edge = 64;
% sig = 10;
% edge = 32;
% sig = 20;
% edge = 48;
step = 1;
w_x=(-edge:step:edge);
w_y=(-edge:step:edge);
[X,Y]=meshgrid(w_x,w_y);
% GaussKernel=exp(-(X.^2+Y.^2)./sig.^2);
ElectrodeMap.GaussKernel = exp(-(X.^2+Y.^2)./sig.^2);

CircMask=zeros(m,n);  
CircMask(circle<=r*r)=1;  %找到圆内的元素，并复制为1
CircMask(circle>r*r)=0;   %找到圆外的元素，并复制为0
CircMaskData.Inside = CircMask;
CircMaskData.Outside = BG_color*double(~CircMaskData.Inside);
CircMaskData.BackGroundColor = uint8(cat(3, BG_color*ones(m,n), BG_color*ones(m,n), BG_color*ones(m,n)));

ElectrodeMap.CircMaskData = CircMaskData;

ElecLocOffSet = 8;
Electrode_x = round([1,2,1,2,1,2,1,2].*repmat((2*r/3),[1,nCh]));
Electrode_x(Electrode_x<r) = Electrode_x(Electrode_x<r)-ElecLocOffSet;
Electrode_x(Electrode_x>r) = Electrode_x(Electrode_x>r)+ElecLocOffSet;
Electrode_y = round([1,1,2,2,3,3,4,4].*repmat((2*r/5),[1,nCh]));
ElectrodeMap.x = Electrode_x;
ElectrodeMap.y = Electrode_y;
ElectrodeMap.ElecLocOffSet = ElecLocOffSet;
ElectrodeMap.ThetaMap = zeros(2*r, 2*r);
ElectrodeMap.AlphaMap = zeros(2*r, 2*r);
ElectrodeMap.BetaMap = zeros(2*r, 2*r);
ElectrodeMap.AllMap = zeros(2*r, 2*r);

ElectrodeMap.ThetaPowerVector = [1,1,1,1,1,1,1,1];
ElectrodeMap.AlphaPowerVector = [1,1,1,1,1,1,1,1];
ElectrodeMap.BetaPowerVector = [1,1,1,1,1,1,1,1];
ElectrodeMap.AllPowerVector = [1,1,1,1,1,1,1,1];
for iCh = 1:nCh
    ElectrodeMap.ThetaMap(ElectrodeMap.y(iCh), ElectrodeMap.x(iCh)) = ElectrodeMap.ThetaPowerVector(iCh);
    ElectrodeMap.AlphaMap(ElectrodeMap.y(iCh), ElectrodeMap.x(iCh)) = ElectrodeMap.AlphaPowerVector(iCh);
    ElectrodeMap.BetaMap(ElectrodeMap.y(iCh), ElectrodeMap.x(iCh)) = ElectrodeMap.BetaPowerVector(iCh);
    ElectrodeMap.AllMap(ElectrodeMap.y(iCh), ElectrodeMap.x(iCh)) = ElectrodeMap.AllPowerVector(iCh);
end
ElectrodeMap.ThetaCovMap = conv2(ElectrodeMap.ThetaMap, ElectrodeMap.GaussKernel, 'same');
ElectrodeMap.AlphaCovMap = conv2(ElectrodeMap.AlphaMap, ElectrodeMap.GaussKernel, 'same');
ElectrodeMap.BetaCovMap = conv2(ElectrodeMap.BetaMap, ElectrodeMap.GaussKernel, 'same');
ElectrodeMap.AllCovMap = conv2(ElectrodeMap.AllMap, ElectrodeMap.GaussKernel, 'same');

% circ_initial = ones(m,n); 
% circ_initial = circ_initial.*255.*CircMaskData.Inside+CircMaskData.Outside;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 设置功率�? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PowerMap.All = 1;
PowerMap.Delta = 0;
PowerMap.Theta = 0;
PowerMap.Alpha = 0;
PowerMap.Beta = 0;
PowerMap.Gamma = 0;
PowerMap.Class = categorical({'Delta','Theta','Alpha','Beta','Gamma'});
PowerMap.Class = categorical(PowerMap.Class,{'Delta','Theta','Alpha','Beta','Gamma'},'Ordinal',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 设置滤波�? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');
[Filter.Coeff_b, Filter.Coeff_a] = tf(Hd);    % Get the transfer function values.
Filter.HistoryOutput = [];   % History y in filter
Filter.TestCoeff_b = Filter.Coeff_b;
Filter.TestCoeff_a = Filter.Coeff_a;
f_Callback = @ButtonCheckFilterResponse_Callback;
f_Callback(handles.ButtonCheckFilterResponse, eventdata, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% 设置EEG Stream矩阵 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DrawMap.Mat = zeros(nCh, SampleRate*DispSec);
DrawMap.CurLoc = 1;
DrawMap.Height = nCh;
DrawMap.Width = SampleRate*DispSec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 设置FFT矩阵 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFTMap.BufferSize = round(SampleRate*IntervalTime.FFTMs/1000);
FFTMap.Buffer = zeros(nCh, FFTMap.BufferSize);
FFTMap.CurLoc = 1;
FFTMap.Height = nCh;
FFTMap.Width = FFTMap.BufferSize;
% fftTest = rand(8,512)
FFTMap.DSR = 2;    % Downsample Rate
FFTMap.Mat = zeros(nCh, floor(SampleRate/FFTMap.DSR/2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 绘图 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.AxesEEGStream, 'YTickLabel', {' ','Fp1','Fp2','F3','F4','P3','P4','O1','O2',' '});

set(handles.AxesEEGStream, 'XLimMode', 'manual');
set(handles.AxesEEGStream, 'YLimMode', 'manual');
set(handles.AxesEEGStream, 'XLim', [1, SampleRate*DispSec]);
set(handles.AxesEEGStream, 'YLim', [0, 90]);

% axes(handles.AxesTopoTheta);
image(handles.AxesTopoTheta, ElectrodeMap.CircMaskData.BackGroundColor); 
hold(handles.AxesTopoTheta, 'on');
imagesc(handles.AxesTopoTheta, ElectrodeMap.ThetaCovMap, 'AlphaData', ElectrodeMap.CircMaskData.Inside);
axis(handles.AxesTopoTheta, 'off');
hold(handles.AxesTopoTheta, 'off');
% axes(handles.AxesTopoAlpha);
image(handles.AxesTopoAlpha, ElectrodeMap.CircMaskData.BackGroundColor); 
hold(handles.AxesTopoAlpha, 'on');
imagesc(handles.AxesTopoAlpha, ElectrodeMap.AlphaCovMap, 'AlphaData', ElectrodeMap.CircMaskData.Inside);
axis(handles.AxesTopoAlpha, 'off');
hold(handles.AxesTopoAlpha, 'off');
% axes(handles.AxesTopoBeta);
image(handles.AxesTopoBeta, ElectrodeMap.CircMaskData.BackGroundColor); 
hold(handles.AxesTopoBeta, 'on');
imagesc(handles.AxesTopoBeta, ElectrodeMap.BetaCovMap, 'AlphaData', ElectrodeMap.CircMaskData.Inside);
axis(handles.AxesTopoBeta, 'off');
hold(handles.AxesTopoBeta, 'off');
% axes(handles.AxesTopoAll);
image(handles.AxesTopoAll, ElectrodeMap.CircMaskData.BackGroundColor); 
hold(handles.AxesTopoAll, 'on');
imagesc(handles.AxesTopoAll, ElectrodeMap.AllCovMap, 'AlphaData', ElectrodeMap.CircMaskData.Inside);
axis(handles.AxesTopoAll, 'off');
hold(handles.AxesTopoAll, 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ZoomCoeff = str2double(get(handles.EditZoomCoeff,'String'));


% Choose default command line output for RTEEGAP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RTEEGAP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RTEEGAP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ButtonOpenSerialPort.
function ButtonOpenSerialPort_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonOpenSerialPort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
instrreset
global ObjSerial;
global SerialPortNum;
global BaudRate;
global SampleRate;
global ZoomCoeff;
global Filter;
global DrawMap;
global BufferSize
global NumTryStartCom
% NumTryStartCom = 0;

% IntervalMs = 200;
% BufferSize = round(33*SampleRate*(IntervalMs/1000));
SerialPortNum = get(handles.EditSerialPort,'String');
ObjSerial=serial(SerialPortNum);                                           %创建串口对象                                                  
ObjSerial.baudrate=BaudRate;                                                 %设置波特�?,缺省9600bit/s
ObjSerial.parity='none';                                                   %无奇偶校�?
ObjSerial.stopbits=1;                                                      %停止�?
ObjSerial.timeout=0.5;                                                     %设置读操作完成时间为1s,缺省10s                                           
ObjSerial.inputbuffersize=BufferSize;                                            %输入缓冲区为32b，缺省�?�为512b


%串口事件回调设置
ObjSerial.BytesAvailableFcnMode='byte';% 设置中断触发事件为�?�bytes-available Event�?
ObjSerial.BytesAvailableFcnCount=BufferSize; % 设置接收缓冲区每收到1800个字节时，触发回调函�?
guidata(hObject,handles);
ObjSerial.BytesAvailableFcn={@mycom, handles}; %得到回调函数句柄
% ObjSerial.BytesAvailableFcn={@instrcallback, ZoomCoeff, Filter, DrawMap}; %得到回调函数句柄
% ObjSerial.BytesAvailableFcnMode='terminator';
% ObjSerial.BytesAvailableFcnCount=10; %输入缓冲区存�?10个字节触发回调函�?
% ObjSerial.BytesAvailableFcn={@EveBytesAvailableFcn,handles};%回调函数的指�?

% handles.TimerHandle_TryStartCom = timer('TimerFcn',{@Timer_TryStartCom,handles},'ExecutionMode', 'fixedDelay', 'Period', 0.5);
    
try
    fopen(ObjSerial);%打开串口
    set(handles.TextCurrentState, 'String', '打开串口成功�?');
catch
    NumTryStartCom = NumTryStartCom+1;
    set(handles.TextCurrentState, 'String', ['正在尝试�?', num2str(NumTryStartCom), '�?']);
%     set(handles.TextCurrentState, 'String', '打开串口失败!');
%     start(handles.TimerHandle_TryStartCom);
end

guidata(hObject,handles);

% --- Executes on button press in ButtonStartRecieveData.
function ButtonStartRecieveData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonStartRecieveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ObjSerial;
% try
%     fwrite(ObjSerial,'b','uchar');                                                               
%     set(handles.TextCurrentState, 'String', '已开始接收数�?!');
% catch
%     set(handles.TextCurrentState, 'String', '发�?�接收指令失�?!');
% end

try
    fwrite(ObjSerial,'b','uchar');                                                                
    set(handles.TextCurrentState, 'String', 'Successful');
catch
    set(handles.TextCurrentState, 'String', 'Failed');
end

guidata(hObject,handles);


% --- Executes on button press in ButtonStopRecieveData.
function ButtonStopRecieveData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonStopRecieveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ObjSerial;

try
    fwrite(ObjSerial,'s','uchar');    
    fclose(ObjSerial);      %关闭串口                                                                
    delete(ObjSerial);      %删除串口对象
%     delete(handles.TimerHandle_TryStartCom)
    set(handles.TextCurrentState, 'String', 'Successful');
catch
    set(handles.TextCurrentState, 'String', 'Failed');
end
guidata(hObject,handles);

function EditSerialPort_Callback(hObject, eventdata, handles)
% hObject    handle to EditSerialPort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSerialPort as text
%        str2double(get(hObject,'String')) returns contents of EditSerialPort as a double
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function EditSerialPort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSerialPort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditZoomCoeff_Callback(hObject, eventdata, handles)
% hObject    handle to EditZoomCoeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function EditZoomCoeff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditZoomCoeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonZoomChange.
function ButtonZoomChange_Callback(hObject, ~, handles)
% hObject    handle to ButtonZoomChange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ZoomCoeff;
ZoomCoeff = str2double(get(handles.EditZoomCoeff,'String'));

function mycom(hObject,eventdata,handles)   %这是串口读取数据后执行的callback函数
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Code By Virter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Filter
global DrawMap
global FFTMap
global ElectrodeMap
global ZoomCoeff
global ObjSerial
global BufferSize;
global SampleRate;
global DispSec;
global SecCnt;
global IntervalTime
global PowerMap
global RecordData
global MotorImagery

if SecCnt<100
    SecCnt = SecCnt+1;
else
    SecCnt = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 常量 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PackageLen = 33;
nCh = 8;
scale_fac_uVolts_per_count = 0.022351744455307063;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 读取数据 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recdta = fread(ObjSerial, BufferSize, 'uchar');   % Receive serial data from OpenBCI (33 Bytes/Package, 256 Package/Second)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 解包 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 滤波 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[UnpackedData, Filter.HistoryOutput] = filter(Filter.Coeff_b, Filter.Coeff_a, UnpackedData, Filter.HistoryOutput, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 记录 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RecordData.RecFlag
    RecordData.NowData = cat(2, RecordData.NowData, UnpackedData);
    RecordData.Time = RecordData.Time + IntervalTime.DispMs/1000;
    RecordData.Sample = RecordData.Sample + size(UnpackedData,2);
    set(handles.TextRecDuration, 'String', [num2str(round(RecordData.Time)), ' s']);

end

if RecordData.EventFlag
    RecordData.EventFlag = 0;
    RecordData.EventTime = [RecordData.EventTime, RecordData.Time];
    RecordData.EventSample = [RecordData.EventSample, RecordData.Sample];
    list=get(handles.ppEventDescribe,'String');
    CurVal=get(handles.ppEventDescribe,'Value');
    RecordData.EventDescribe = [RecordData.EventDescribe, list{CurVal}];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 运动想象 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if MotorImagery.RecordFlag

    MotorImagery.RecordData = cat(2, MotorImagery.RecordData, UnpackedData);
    MotorImagery.RecordRawData = cat(2, MotorImagery.RecordRawData, UnpackedDataRaw);
    MotorImagery.Sample.Cnt = MotorImagery.Sample.Cnt + size(UnpackedData,2);
    MotorImagery.EpochCnt = MotorImagery.EpochCnt + 1;
    if (MotorImagery.EpochCnt >= MotorImagery.PrepareTime/(IntervalTime.DispMs/1000)) && (MotorImagery.StageFlag == 0)
        MotorImagery.EpochCnt = 0;
        MotorImagery.StageFlag = 1;
        MotorImagery.Sample.Stage1 = MotorImagery.Sample.Cnt;
        MotorImagery.State = 'Preparing';
        %%%%% 显示图像 %%%%%
        if logical(get(handles.CheckBoxEnablePredict, 'Value')) && ~logical(get(handles.CheckBoxEnableCue, 'Value'))
            if get(handles.CheckBoxUseBigCue, 'Value')
                image(handles.AxesImageBigCue, MotorImagery.NoCueImage);
                axis(handles.AxesImageBigCue, 'off');
            end
            image(handles.AxesImageSmallCue, MotorImagery.NoCueImage);
            axis(handles.AxesImageSmallCue, 'off');              
        else
            CurrentClass = MotorImagery.ClassOrder(MotorImagery.CurrentCnt+1);  % Cnt+1
            if get(handles.CheckBoxUseBigCue, 'Value')
                image(handles.AxesImageBigCue, MotorImagery.ClassImage{CurrentClass});
                axis(handles.AxesImageBigCue, 'off');
            end
            image(handles.AxesImageSmallCue, MotorImagery.ClassImage{CurrentClass});
            axis(handles.AxesImageSmallCue, 'off');
        end
        %%%%%%%%%%%%%%%%%%%%        
    end
    if (MotorImagery.EpochCnt >= MotorImagery.CueTime/(IntervalTime.DispMs/1000)) && (MotorImagery.StageFlag == 1)
        MotorImagery.EpochCnt = 0;
        MotorImagery.StageFlag = 2;
        MotorImagery.Sample.Stage2 = MotorImagery.Sample.Cnt;
        MotorImagery.State = 'Imaging';
        %%%%% 显示空白 %%%%%
        if get(handles.CheckBoxUseBigCue, 'Value')
            image(handles.AxesImageBigCue, cat(3,ones(5,5),ones(5,5),ones(5,5)));
            axis(handles.AxesImageBigCue, 'off');
        end
        image(handles.AxesImageSmallCue, cat(3,ones(5,5),ones(5,5),ones(5,5)));
        axis(handles.AxesImageSmallCue, 'off');        
        %%%%%%%%%%%%%%%%%%%%        
    end
    if (MotorImagery.EpochCnt >= MotorImagery.ImageLastTime/(IntervalTime.DispMs/1000)) && (MotorImagery.StageFlag == 2)
        MotorImagery.EpochCnt = 0;
        MotorImagery.StageFlag = 3;
        MotorImagery.Sample.Stage3 = MotorImagery.Sample.Cnt;
        MotorImagery.State = 'Stopping';
        %%%%% 显示休息 %%%%%
        if logical(get(handles.CheckBoxEnablePredict, 'Value'))
            ImageData(1,:,:) = MotorImagery.RecordData(:, MotorImagery.Sample.Stage2+1:MotorImagery.Sample.Stage2+1000);
            SelFlag = get(handles.CheckBoxEnableFeaSel, 'Value');
            PredictLabel = PredictSingleTrail(MotorImagery.ParaImagery, ImageData, SelFlag, MotorImagery.ParaImagery.CSP_Config);
            PredictLabel = double(PredictLabel);
            image(handles.AxesImageBigCue, MotorImagery.ClassImage{PredictLabel});
            axis(handles.AxesImageBigCue, 'off');
            image(handles.AxesImageSmallCue, MotorImagery.ClassImage{PredictLabel});
            axis(handles.AxesImageSmallCue, 'off');            
        else    
            if get(handles.CheckBoxUseBigCue, 'Value')
                image(handles.AxesImageBigCue, cat(3,zeros(5,5),zeros(5,5),zeros(5,5)));
                axis(handles.AxesImageBigCue, 'off');
            end
            image(handles.AxesImageSmallCue, cat(3,zeros(5,5),zeros(5,5),zeros(5,5)));
            axis(handles.AxesImageSmallCue, 'off');   
        end
        %%%%%%%%%%%%%%%%%%%%        
    end
    if (MotorImagery.EpochCnt >= MotorImagery.AllTime/(IntervalTime.DispMs/1000)) && (MotorImagery.StageFlag == 3)
        MotorImagery.EpochCnt = 0;
        MotorImagery.StageFlag = 4;
        MotorImagery.RecordFlag = 0;
        MotorImagery.Sample.Stage4 = MotorImagery.Sample.Cnt;
        MotorImagery.State = 'Stopping';
        %%%%% 保存数据 %%%%%        
        DrawMapTrans = MotorImagery.RecordData./ZoomCoeff;
        interval = 10;
        ch_offset = (interval:interval:interval*nCh)';
        ch_offset = repmat(ch_offset, [1,size(DrawMapTrans,2)]);
        DrawMapTrans = DrawMapTrans+ch_offset;
        DrawMapTrans(DrawMapTrans>90) = 90;
        DrawMapTrans(DrawMapTrans<0) = 0;
        
        plot(handles.AxesEEGStream, DrawMapTrans');
        set(handles.AxesEEGStream, 'YTickLabel', {' ','Fp1','Fp2','F3','F4','P3','P4','O1','O2',' '});

        if get(handles.CheckBoxUseBigCue, 'Value')
            plot(handles.AxesImageBigCue, 0,0);
            set(handles.AxesImageBigCue, 'Visible', 'off');
            set(handles.uipanel1, 'Visible', 'on');
            set(handles.uipanel2, 'Visible', 'on');
            set(handles.uipanel3, 'Visible', 'on');
        end
        
        answer = questdlg('Do you want to save this trial?', ...
                          'Option', ...
                          'Save','Do Not Save','Save');
        switch answer
            case 'Save'
                MotorImagery.RemainEpoch = MotorImagery.RemainEpoch - 1;
                set(handles.TextRemainEpoch, 'String', num2str(MotorImagery.RemainEpoch));
                MotorImagery.CurrentCnt = MotorImagery.CurrentCnt + 1;
                save([MotorImagery.ImageFilePath, num2str(MotorImagery.CurrentCnt), '.mat'], 'MotorImagery');
            case 'Do Not Save'

        end
        
        %%%%%%%%%%%%%%%%%%%%
    end
    
    if MotorImagery.StageFlag > 3
        MotorImagery.EpochCnt = 0;
        MotorImagery.RecordFlag = 0;
        MotorImagery.StageFlag =  0;
        MotorImagery.Sample.Cnt = 0;                                                                    % 样本点数
        MotorImagery.Sample.Stage1 = 0;
        MotorImagery.Sample.Stage2 = 0;
        MotorImagery.Sample.Stage3 = 0;
        MotorImagery.Sample.Stage4 = 0;
        MotorImagery.RecordData = [];
        MotorImagery.RecordRawData = [];
        MotorImagery.State = 'Not Start';
    end
    
    set(handles.TextImageExperimentState, 'String', MotorImagery.State);    % 用此句实时更新状�?
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 绘制 EEG Stream %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DrawMap.CurLoc > DrawMap.Width 
    DrawMap.CurLoc = 1; 
end

DataLength = size(UnpackedData,2);
if DataLength+DrawMap.CurLoc > DrawMap.Width
    dataPart = (DrawMap.Width-DrawMap.CurLoc)+1;
    DrawMap.Mat(:, DrawMap.CurLoc:end) = UnpackedData(:, 1:dataPart);
    DrawMap.Mat(:,1:(DataLength-dataPart)) = UnpackedData(:, dataPart+1:end);
    DrawMap.CurLoc = (DataLength-dataPart)+1;
else
    DrawMap.Mat(:, DrawMap.CurLoc:DrawMap.CurLoc+DataLength-1) = UnpackedData;
    DrawMap.CurLoc = DrawMap.CurLoc + DataLength;
end

if DrawMap.CurLoc > DrawMap.Width 
    DrawMap.CurLoc = 1; 
end

DrawMapTrans = DrawMap.Mat./ZoomCoeff;
interval = 10;
ch_offset = (interval:interval:interval*nCh)';
ch_offset = repmat(ch_offset, [1,size(DrawMapTrans,2)]);
DrawMapTrans = DrawMapTrans+ch_offset;

DrawMapTrans(DrawMapTrans>90) = 90;
DrawMapTrans(DrawMapTrans<0) = 0;

set(handles.TextUnpackedData, 'String', num2str(DataLength));
% save(['I:\STUDY\MATLAB上位机\2\tmpSerialData\SerialData_', num2str(SecCnt), '.mat'], 'recdta', 'SegIdx', 'PackArr', 'UnpackedData', 'Filter', 'DrawMap');

plot(handles.AxesEEGStream, DrawMapTrans');
hold(handles.AxesEEGStream, 'on');
plot(handles.AxesEEGStream,[DrawMap.CurLoc,DrawMap.CurLoc],[0,90],'r');
hold(handles.AxesEEGStream, 'off');
set(handles.AxesEEGStream, 'YTickLabel', {' ','Fp1','Fp2','F3','F4','P3','P4','O1','O2',' '});

% set(handles.AxesEEGStream, 'XLimMode', 'manual');
% set(handles.AxesEEGStream, 'YLimMode', 'manual');
set(handles.AxesEEGStream, 'XLim', [1, SampleRate*DispSec]);
set(handles.AxesEEGStream, 'YLim', [0, 90]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DataLength > FFTMap.Width
    FFTMap.Buffer = UnpackedData(:, end-FFTMap.Width+1:end);
else
    FFTMap.Buffer(:, 1:(end-DataLength)) = FFTMap.Buffer(:, DataLength+1:end);
    FFTMap.Buffer(:, (end-DataLength+1):end) = UnpackedData;
end

FFTBufferDS = FFTMap.Buffer(:,1:FFTMap.DSR:end);   % Downsampled by 2 times
[x,fft_mold] = fftmold(FFTBufferDS,SampleRate/FFTMap.DSR,SampleRate/FFTMap.DSR);
plot(handles.AxesFFT, x(1:40), fft_mold(:,1:40));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 绘制功率�? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PowerMap.All = sum(sum(fft_mold(:,1:45).^2));
PowerMap.Delta = sum(sum(fft_mold(:,1:3).^2))/PowerMap.All;
PowerMap.Theta = sum(sum(fft_mold(:,4:7).^2))/PowerMap.All;
PowerMap.Alpha = sum(sum(fft_mold(:,8:13).^2))/PowerMap.All;
PowerMap.Beta = sum(sum(fft_mold(:,14:30).^2))/PowerMap.All;
PowerMap.Gamma = sum(sum(fft_mold(:,31:45).^2))/PowerMap.All;
bar(handles.AxesPower, PowerMap.Class, [PowerMap.Delta, PowerMap.Theta, PowerMap.Alpha, PowerMap.Beta, PowerMap.Gamma]);
set(handles.AxesPower, 'YLim', [0,1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 脑电地形�? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ElectrodeMap.ThetaPowerVector = sum(fft_mold(:,4:7).^2,2)';
ElectrodeMap.AlphaPowerVector = sum(fft_mold(:,8:13).^2,2)';
ElectrodeMap.BetaPowerVector = sum(fft_mold(:,14:30).^2,2)';
ElectrodeMap.AllPowerVector = sum(fft_mold(:,4:30).^2,2)';
for iCh = 1:nCh
    ElectrodeMap.ThetaMap(ElectrodeMap.y(iCh), ElectrodeMap.x(iCh)) = ElectrodeMap.ThetaPowerVector(iCh);
    ElectrodeMap.AlphaMap(ElectrodeMap.y(iCh), ElectrodeMap.x(iCh)) = ElectrodeMap.AlphaPowerVector(iCh);
    ElectrodeMap.BetaMap(ElectrodeMap.y(iCh), ElectrodeMap.x(iCh)) = ElectrodeMap.BetaPowerVector(iCh);
    ElectrodeMap.AllMap(ElectrodeMap.y(iCh), ElectrodeMap.x(iCh)) = ElectrodeMap.AllPowerVector(iCh);
end
ElectrodeMap.ThetaCovMap = conv2(ElectrodeMap.ThetaMap, ElectrodeMap.GaussKernel, 'same');
ElectrodeMap.AlphaCovMap = conv2(ElectrodeMap.AlphaMap, ElectrodeMap.GaussKernel, 'same');
ElectrodeMap.BetaCovMap = conv2(ElectrodeMap.BetaMap, ElectrodeMap.GaussKernel, 'same');
ElectrodeMap.AllCovMap = conv2(ElectrodeMap.AllMap, ElectrodeMap.GaussKernel, 'same');

% axes(handles.AxesTopoTheta);
image(handles.AxesTopoTheta, ElectrodeMap.CircMaskData.BackGroundColor); 
hold(handles.AxesTopoTheta, 'on');
imagesc(handles.AxesTopoTheta, ElectrodeMap.ThetaCovMap, 'AlphaData', ElectrodeMap.CircMaskData.Inside);
axis(handles.AxesTopoTheta, 'off');
hold(handles.AxesTopoTheta, 'off');
% axes(handles.AxesTopoAlpha);
image(handles.AxesTopoAlpha, ElectrodeMap.CircMaskData.BackGroundColor); 
hold(handles.AxesTopoAlpha, 'on');
imagesc(handles.AxesTopoAlpha, ElectrodeMap.AlphaCovMap, 'AlphaData', ElectrodeMap.CircMaskData.Inside);
axis(handles.AxesTopoAlpha, 'off');
hold(handles.AxesTopoAlpha, 'off');
% axes(handles.AxesTopoBeta);
image(handles.AxesTopoBeta, ElectrodeMap.CircMaskData.BackGroundColor); 
hold(handles.AxesTopoBeta, 'on');
imagesc(handles.AxesTopoBeta, ElectrodeMap.BetaCovMap, 'AlphaData', ElectrodeMap.CircMaskData.Inside);
axis(handles.AxesTopoBeta, 'off');
hold(handles.AxesTopoBeta, 'off');
% axes(handles.AxesTopoAll);
image(handles.AxesTopoAll, ElectrodeMap.CircMaskData.BackGroundColor); 
hold(handles.AxesTopoAll, 'on');
imagesc(handles.AxesTopoAll, ElectrodeMap.AllCovMap, 'AlphaData', ElectrodeMap.CircMaskData.Inside);
axis(handles.AxesTopoAll, 'off');
hold(handles.AxesTopoAll, 'off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(handles.TextLoad, 'String', [num2str(double(vpa(toc*1000/IntervalTime.DispMs,4))*100), '%']);
set(handles.TextLoad, 'String', [num2str(round(toc*1000/IntervalTime.DispMs*100)), '%']);

% --- Executes during object creation, after setting all properties.
function AxesTopoTheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AxesTopoTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate AxesTopoTheta
% set(hObject, 'XColor', [240/255,240/255,240/255]);
% set(hObject, 'YColor', [240/255,240/255,240/255]);
% axis(hObject, 'off');



function EditOrder_Callback(hObject, eventdata, handles)
% hObject    handle to EditOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditOrder as text
%        str2double(get(hObject,'String')) returns contents of EditOrder as a double


% --- Executes during object creation, after setting all properties.
function EditOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonChangeFilter.
function ButtonChangeFilter_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonChangeFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Filter
Filter.Coeff_a = Filter.TestCoeff_a;
Filter.Coeff_b = Filter.TestCoeff_b;
Filter.HistoryOutput = [];   % History y in filter

function EditFc1_Callback(hObject, eventdata, handles)
% hObject    handle to EditFc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditFc1 as text
%        str2double(get(hObject,'String')) returns contents of EditFc1 as a double


% --- Executes during object creation, after setting all properties.
function EditFc1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditFc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditFc2_Callback(hObject, eventdata, handles)
% hObject    handle to EditFc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditFc2 as text
%        str2double(get(hObject,'String')) returns contents of EditFc2 as a double


% --- Executes during object creation, after setting all properties.
function EditFc2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditFc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonCheckFilterResponse.
function ButtonCheckFilterResponse_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonCheckFilterResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SampleRate
global Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filter Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = str2double(get(handles.EditOrder,'String'));
Fc1 = str2double(get(handles.EditFc1,'String'));
Fc2 = str2double(get(handles.EditFc2,'String'));
Fs = SampleRate;
% FilterType = get(handles.ppFilterType,'String');
list=get(handles.ppFilterType,'String');
val1=get(handles.ppFilterType,'Value');
FilterType = list{val1};

if strcmp(FilterType, 'IIR')  
    h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
    Hd = design(h, 'butter');
    [Filter.TestCoeff_b, Filter.TestCoeff_a] = tf(Hd);    % Get the transfer function values.
    Response = freqz(Filter.TestCoeff_b,Filter.TestCoeff_a,round(Fs/2));
else
    flag = 'scale';  % Sampling Flag
    win = hann(N+1);
    Filter.TestCoeff_b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
    Filter.TestCoeff_a = 1;
    Response = freqz(Filter.TestCoeff_b,1,round(Fs/2));
end
% Filter.HistoryOutput = [];   % History y in filter
hAx = plotyy(handles.AxesFilterResponse,1:round(Fs/2),20*log10(abs(Response)),1:round(Fs/2),unwrap(angle(Response))*180/pi);
xlabel(handles.AxesFilterResponse, 'Frequency (Hz)');
ylabel(hAx(1), 'Magnitude (dB)');
ylabel(hAx(2), 'Angle (degrees)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppFilterType.
function ppFilterType_Callback(hObject, eventdata, handles)
% hObject    handle to ppFilterType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppFilterType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppFilterType


% --- Executes during object creation, after setting all properties.
function ppFilterType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppFilterType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSavePath_Callback(hObject, eventdata, handles)
% hObject    handle to EditSavePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSavePath as text
%        str2double(get(hObject,'String')) returns contents of EditSavePath as a double


% --- Executes during object creation, after setting all properties.
function EditSavePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSavePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonRecordData.
function ButtonRecordData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRecordData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RecordData
if RecordData.RecFlag
    RecordData.RecFlag = 0;
    set(handles.ButtonRecordData, 'String', 'Start Record');
    set(handles.TextRecDuration, 'String', 'Saving......');
    RecordData.Path = get(handles.EditSavePath,'String');
    save(RecordData.Path, 'RecordData');
    set(handles.TextRecDuration, 'String', 'Finished!');
    RecordData.NowData = [];
else
    RecordData.RecFlag = 1;
    RecordData.NowData = [];
    RecordData.EventFlag = 0;
    RecordData.Path = get(handles.EditSavePath,'String');
    RecordData.Time = 0;
    RecordData.Sample = 0;
    RecordData.EventTime = [];
    RecordData.EventSample = [];
    RecordData.EventDescribe = {};
    set(handles.ButtonRecordData, 'String', 'Stop Record');
    set(handles.TextRecDuration, 'String', '0 s');
end

% --- Executes on selection change in ppEventDescribe.
function ppEventDescribe_Callback(hObject, eventdata, handles)
% hObject    handle to ppEventDescribe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppEventDescribe contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppEventDescribe


% --- Executes during object creation, after setting all properties.
function ppEventDescribe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppEventDescribe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonAddEvent.
function ButtonAddEvent_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonAddEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RecordData
RecordData.EventFlag = 1;


% --- Executes on button press in ButtonAddImageEvent.
function ButtonAddImageEvent_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonAddImageEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MotorImagery
ListContent = get(handles.ListBoxImageEvent, 'String');
ListCurVal = get(handles.ListBoxImageEvent, 'Value');
if ListCurVal<numel(ListContent)
    ListContent(ListCurVal+2:numel(ListContent)+1) = ListContent(ListCurVal+1:numel(ListContent));
    ListContent{ListCurVal+1} = get(handles.EditImageEvent, 'String');
    set(handles.ListBoxImageEvent, 'String', ListContent);
else
    ListContent{ListCurVal+1} = get(handles.EditImageEvent, 'String');
    set(handles.ListBoxImageEvent, 'String', ListContent);
end

MotorImagery.ClassName = get(handles.ListBoxImageEvent, 'String');
MotorImagery.ClassEpoch = str2double(get(handles.EditClassEpoch, 'String'));
MotorImagery.numClass = numel(MotorImagery.ClassName);
MotorImagery.AllEpoch = MotorImagery.numClass*MotorImagery.ClassEpoch;                % 总试验次�?
MotorImagery.RemainEpoch = MotorImagery.AllEpoch;                                           % 剩余次数
set(handles.TextRemainEpoch, 'String', num2str(MotorImagery.RemainEpoch));
ClassOrder = repmat(1:MotorImagery.numClass,[1 MotorImagery.ClassEpoch]);
ClassOrder = ClassOrder(randperm(numel(ClassOrder)));
MotorImagery.ClassOrder = ClassOrder;

function EditImageEvent_Callback(hObject, eventdata, handles)
% hObject    handle to EditImageEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditImageEvent as text
%        str2double(get(hObject,'String')) returns contents of EditImageEvent as a double


% --- Executes during object creation, after setting all properties.
function EditImageEvent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditImageEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ListBoxImageEvent.
function ListBoxImageEvent_Callback(hObject, eventdata, handles)
% hObject    handle to ListBoxImageEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListBoxImageEvent contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListBoxImageEvent


% --- Executes during object creation, after setting all properties.
function ListBoxImageEvent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListBoxImageEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonDeleteImageEvent.
function ButtonDeleteImageEvent_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonDeleteImageEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MotorImagery
ListContent = get(handles.ListBoxImageEvent, 'String');
ListCurVal = get(handles.ListBoxImageEvent, 'Value');
if numel(ListContent)<=1
    warndlg('The number of the event must Ggreater than 1!','WARNING','modal');
else
    ListContent(ListCurVal) = [];
    set(handles.ListBoxImageEvent, 'Value', 1);
    set(handles.ListBoxImageEvent, 'String', ListContent);
end
MotorImagery.ClassName = get(handles.ListBoxImageEvent, 'String');
MotorImagery.ClassEpoch = str2double(get(handles.EditClassEpoch, 'String'));
MotorImagery.numClass = numel(MotorImagery.ClassName);
MotorImagery.AllEpoch = MotorImagery.numClass*MotorImagery.ClassEpoch;                % 总试验次�?
MotorImagery.RemainEpoch = MotorImagery.AllEpoch;                                           % 剩余次数
set(handles.TextRemainEpoch, 'String', num2str(MotorImagery.RemainEpoch));
ClassOrder = repmat(1:MotorImagery.numClass,[1 MotorImagery.ClassEpoch]);
ClassOrder = ClassOrder(randperm(numel(ClassOrder)));
MotorImagery.ClassOrder = ClassOrder;

function EditPrepareTime_Callback(hObject, eventdata, handles)
% hObject    handle to EditPrepareTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditPrepareTime as text
%        str2double(get(hObject,'String')) returns contents of EditPrepareTime as a double


% --- Executes during object creation, after setting all properties.
function EditPrepareTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditPrepareTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditImageLastTime_Callback(hObject, eventdata, handles)
% hObject    handle to EditImageLastTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditImageLastTime as text
%        str2double(get(hObject,'String')) returns contents of EditImageLastTime as a double


% --- Executes during object creation, after setting all properties.
function EditImageLastTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditImageLastTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditAllTime_Callback(hObject, eventdata, handles)
% hObject    handle to EditAllTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditAllTime as text
%        str2double(get(hObject,'String')) returns contents of EditAllTime as a double


% --- Executes during object creation, after setting all properties.
function EditAllTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditAllTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditImageFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to EditImageFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditImageFilePath as text
%        str2double(get(hObject,'String')) returns contents of EditImageFilePath as a double


% --- Executes during object creation, after setting all properties.
function EditImageFilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditImageFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditClassEpoch_Callback(hObject, eventdata, handles)
% hObject    handle to EditClassEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditClassEpoch as text
%        str2double(get(hObject,'String')) returns contents of EditClassEpoch as a double


% --- Executes during object creation, after setting all properties.
function EditClassEpoch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditClassEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonStartExperiment.
function ButtonStartExperiment_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonStartExperiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MotorImagery

% MotorImagery.RemainEpoch = MotorImagery.RemainEpoch - 1;
if MotorImagery.RemainEpoch-1 < 0
    warndlg('This experiment is finished!','WARNING','modal');
else
    MotorImagery.PrepareTime = str2double(get(handles.EditPrepareTime, 'String'));              % 准备时间
    MotorImagery.ImageLastTime = str2double(get(handles.EditImageLastTime, 'String'));          % 持续时间
    MotorImagery.CueTime = str2double(get(handles.EditCueTime, 'String'));                      % 提示时间
    MotorImagery.AllTime = str2double(get(handles.EditAllTime, 'String'));                      % 末尾时间
    MotorImagery.ClassEpoch = str2double(get(handles.EditClassEpoch, 'String'));
    MotorImagery.RemainEpoch = str2double(get(handles.TextRemainEpoch, 'String'));
    MotorImagery.ImageFilePath = get(handles.EditImageFilePath, 'String');                      % 运动想象数据保存路径
    save([MotorImagery.CurrentPath, '\Data\Configure\Configure_AutoSave.mat'], 'MotorImagery');
%     set(handles.TextRemainEpoch, 'String', num2str(MotorImagery.RemainEpoch));
    MotorImagery.RecordFlag = 1;
    MotorImagery.EpochCnt = 0;
%     MotorImagery.CurrentCnt = MotorImagery.CurrentCnt + 1;
    if get(handles.CheckBoxUseBigCue, 'Value')
        set(handles.AxesImageBigCue, 'Visible', 'on');
        image(handles.AxesImageBigCue, MotorImagery.PrepareImage);
        axis(handles.AxesImageBigCue, 'off');
        set(handles.uipanel1, 'Visible', 'off');
        set(handles.uipanel2, 'Visible', 'off');
        set(handles.uipanel3, 'Visible', 'off');
    end
    image(handles.AxesImageSmallCue, MotorImagery.PrepareImage);
    axis(handles.AxesImageSmallCue, 'off');
    sound(MotorImagery.ImageSound.SoundCue, MotorImagery.ImageSound.SoundFs);   % 播放音乐
end

function EditImageConfigureFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to EditImageConfigureFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditImageConfigureFilePath as text
%        str2double(get(hObject,'String')) returns contents of EditImageConfigureFilePath as a double


% --- Executes during object creation, after setting all properties.
function EditImageConfigureFilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditImageConfigureFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonSaveImageConfig.
function ButtonSaveImageConfig_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonSaveImageConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MotorImagery
MotorImagery.PrepareTime = str2double(get(handles.EditPrepareTime, 'String'));              % 准备时间
MotorImagery.ImageLastTime = str2double(get(handles.EditImageLastTime, 'String'));          % 持续时间
MotorImagery.CueTime = str2double(get(handles.EditCueTime, 'String'));                      % 提示时间
MotorImagery.AllTime = str2double(get(handles.EditAllTime, 'String'));                      % 末尾时间
MotorImagery.ClassEpoch = str2double(get(handles.EditClassEpoch, 'String'));
MotorImagery.RemainEpoch = str2double(get(handles.TextRemainEpoch, 'String'));
MotorImagery.ImageFilePath = get(handles.EditImageFilePath, 'String');                      % 运动想象数据保存路径

save(get(handles.EditImageConfigureFilePath, 'String'), 'MotorImagery');


% --- Executes on button press in ButtonLoadImageConfig.
function ButtonLoadImageConfig_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoadImageConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MotorImagery
load(get(handles.EditImageConfigureFilePath, 'String'), 'MotorImagery');

set(handles.EditPrepareTime, 'String', num2str(MotorImagery.PrepareTime));
set(handles.EditImageLastTime, 'String', num2str(MotorImagery.ImageLastTime));
set(handles.EditCueTime,  'String', num2str(MotorImagery.CueTime));
set(handles.EditAllTime, 'String', num2str(MotorImagery.AllTime));
set(handles.EditImageFilePath, 'String', MotorImagery.ImageFilePath);
% set(handles.EditImageConfigureFilePath, 'String', MotorImagery.ImageConfigureFilePath);
set(handles.ListBoxImageEvent, 'String', MotorImagery.ClassName);
set(handles.EditClassEpoch, 'String', num2str(MotorImagery.ClassEpoch));
set(handles.TextImageExperimentState, 'String', MotorImagery.State);
set(handles.TextRemainEpoch, 'String', num2str(MotorImagery.RemainEpoch));
set(handles.CheckBoxUseBigCue, 'Value', MotorImagery.UseBigCue);

% --- Executes on button press in CheckBoxUseBigCue.
function CheckBoxUseBigCue_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxUseBigCue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxUseBigCue


% --- Executes on button press in CheckBoxEnableCue.
function CheckBoxEnableCue_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxEnableCue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxEnableCue



function EditCueTime_Callback(hObject, eventdata, handles)
% hObject    handle to EditCueTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditCueTime as text
%        str2double(get(hObject,'String')) returns contents of EditCueTime as a double


% --- Executes during object creation, after setting all properties.
function EditCueTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditCueTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CheckBoxEnablePredict.
function CheckBoxEnablePredict_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxEnablePredict (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxEnablePredict

if get(hObject,'Value')
    set(handles.CheckBoxUseBigCue, 'Value', 1);
    set(handles.CheckBoxUseBigCue, 'Enable', 'off');
    set(handles.CheckBoxEnableFeaSel, 'Enable', 'on');
    set(handles.CheckBoxEnableCue, 'Enable', 'on');
else
    set(handles.CheckBoxUseBigCue, 'Enable', 'on');
    set(handles.CheckBoxEnableFeaSel, 'Enable', 'off');
    set(handles.CheckBoxEnableCue, 'Enable', 'off');
end

% --- Executes on button press in CheckBoxEnableFeaSel.
function CheckBoxEnableFeaSel_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxEnableFeaSel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxEnableFeaSel



function EditPredictFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to EditPredictFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditPredictFilePath as text
%        str2double(get(hObject,'String')) returns contents of EditPredictFilePath as a double


% --- Executes during object creation, after setting all properties.
function EditPredictFilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditPredictFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonLoadPredictFile.
function ButtonLoadPredictFile_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoadPredictFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MotorImagery
[filename,filepath]=uigetfile('*.mat','Please select the file');
load([filepath,filename], 'ParaImagery');
MotorImagery.ParaImagery = ParaImagery;
set(handles.CheckBoxEnablePredict, 'Enable', 'on');
% set(handles.CheckBoxEnableFeaSel, 'Enable', 'on');
% set(handles.CheckBoxEnableCue, 'Enable', 'on');
