function varargout = MTP(varargin)
% MTP MATLAB code for MTP.fig
%      MTP, by itself, creates a new MTP or raises the existing
%      singleton*.
%
%      H = MTP returns the handle to a new MTP or the handle to
%      the existing singleton*.
%
%      MTP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MTP.M with the given input arguments.
%
%      MTP('Property','Value',...) creates a new MTP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MTP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MTP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MTP

% Last Modified by GUIDE v2.5 09-May-2022 13:03:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MTP_OpeningFcn, ...
                   'gui_OutputFcn',  @MTP_OutputFcn, ...
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


% --- Executes just before MTP is made visible.
function MTP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MTP (see VARARGIN)

% Choose default command line output for MTP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global TrainSetting


CurrentPath = cd;
ThisToolRootPath = [CurrentPath, '\'];
addpath(genpath(ThisToolRootPath));

TrainSetting.TrainingDataFold = ThisToolRootPath;
TrainSetting.TrainedModelSavePath = [ThisToolRootPath, 'Data\Configure\TrainedModel.mat'];
TrainSetting.ImageData = [];
TrainSetting.ImageLabel = [];

set(handles.EditTrainParaFilePath, 'String', TrainSetting.TrainedModelSavePath);

% UIWAIT makes MTP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MTP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EditTrainDataPath_Callback(hObject, eventdata, handles)
% hObject    handle to EditTrainDataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditTrainDataPath as text
%        str2double(get(hObject,'String')) returns contents of EditTrainDataPath as a double


% --- Executes during object creation, after setting all properties.
function EditTrainDataPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditTrainDataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CheckBoxFeaSel.
function CheckBoxFeaSel_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxFeaSel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxFeaSel


% --- Executes on button press in CheckBoxEnableParallel.
function CheckBoxEnableParallel_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxEnableParallel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxEnableParallel


% --- Executes on selection change in ppCSPMode.
function ppCSPMode_Callback(hObject, eventdata, handles)
% hObject    handle to ppCSPMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppCSPMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppCSPMode


% --- Executes during object creation, after setting all properties.
function ppCSPMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppCSPMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditNumFold_Callback(hObject, eventdata, handles)
% hObject    handle to EditNumFold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditNumFold as text
%        str2double(get(hObject,'String')) returns contents of EditNumFold as a double


% --- Executes during object creation, after setting all properties.
function EditNumFold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditNumFold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonStartTrainingAllData.
function ButtonStartTrainingAllData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonStartTrainingAllData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global TrainSetting

SelFlag = get(handles.CheckBoxFeaSel, 'Value');    % 是否进行数据特征选择，0：无特征选择，1：特征选择
ParallelFlag = get(handles.CheckBoxEnableParallel, 'Value');   % 是否启用并行，0：不启用并行，1：启用并行
list=get(handles.ppCSPMode,'String');
CurVal=get(handles.ppCSPMode,'Value');
switch list{CurVal}
    case 'CSP'
        CSP_Config.Mode = 2;    % CSP模式，1：divCSP, 2：CSP
    case 'div-CSP'
        CSP_Config.Mode = 1;    % CSP模式，1：divCSP, 2：CSP
    otherwise
        error('No Such Option');
end

% CSP_Config.Mode = 2;    % CSP模式，1：divCSP, 2：CSP
CSP_Config.Wcsp = [];
nClass = 2;     % 分类数目

if ~isempty(TrainSetting.ImageData)
    ParaImagery = GetParaImagery(TrainSetting.ImageData, TrainSetting.ImageLabel, SelFlag, ParallelFlag, CSP_Config, nClass);
    ParaImagery.CSP_Config = CSP_Config;
else
    error('Please load the training data');
end

SavePath = get(handles.EditTrainParaFilePath, 'String');
save(SavePath, 'ParaImagery');
set(handles.TextState, 'String', ['Training is completed! Model has been saved!']);

% --- Executes on button press in ButtonStartCrossValidation.
function ButtonStartCrossValidation_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonStartCrossValidation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global TrainSetting

SelFlag = get(handles.CheckBoxFeaSel, 'Value');    % 是否进行数据特征选择，0：无特征选择，1：特征选择
ParallelFlag = get(handles.CheckBoxEnableParallel, 'Value');   % 是否启用并行，0：不启用并行，1：启用并行
list=get(handles.ppCSPMode,'String');
CurVal=get(handles.ppCSPMode,'Value');
switch list{CurVal}
    case 'CSP'
        CSP_Config.Mode = 2;    % CSP模式，1：divCSP, 2：CSP
    case 'div-CSP'
        CSP_Config.Mode = 1;    % CSP模式，1：divCSP, 2：CSP
    otherwise
        error('No Such Options');
end

% CSP_Config.Mode = 2;    % CSP模式，1：divCSP, 2：CSP
CSP_Config.Wcsp = [];
nFold = str2double(get(handles.EditNumFold, 'String'));      % 交叉验证折数
nClass = 2;     % 分类数目

accMat = [];
indices = crossvalind('Kfold',numel(find(double(TrainSetting.ImageLabel)==1)),nFold);
% load divCSP_Selection_5Folds.mat indices

for iFold = 1:nFold
    tic;
    n_test = (indices == iFold);
    n_train = ~n_test;
    [Data.Train, Data.Test, Data.TrainLabel, Data.TestLabel] = GetFoldData(TrainSetting.ImageData, TrainSetting.ImageLabel, n_train, n_test, nClass);
        
    if SelFlag
        indices_train = crossvalind('Kfold',numel(find(double(Data.TrainLabel)==1)),2);
        [DataSel.Train, DataSel.Test, DataSel.TrainLabel, DataSel.TestLabel] = GetFoldData(Data.Train, Data.TrainLabel, ~(indices_train == 1), (indices_train == 1), nClass);
        [~, ~, ~, ~, Wcsp_Sel, PatFeature_Sel] = FBNN_Excute_CPIII(DataSel.Train,DataSel.TrainLabel,DataSel.Test,DataSel.TestLabel,ParallelFlag,CSP_Config,1);
        [accMat{1}, kappaMat{1}] = FBNN_Excute_FeatureSelection(PatFeature_Sel,ParallelFlag,1);
        TFMat = PlotTFMap_All_Return(accMat);

        s = pcolor(handles.AxesFeaSelVisualize, 0.5:0.5:4, 6:2:34, TFMat{1}(2:end-1, :));
        s.FaceColor = 'interp';
        s.LineStyle = 'none';
        set(handles.AxesFeaSelVisualize,'YDir','normal'); 
        xlabel(handles.AxesFeaSelVisualize, 'Time(s)');
        ylabel(handles.AxesFeaSelVisualize, 'Frequency(Hz)');
    end
    [~, ~, ~, ~, Wcsp, PatFeature] = FBNN_Excute_CPIII(Data.Train,Data.TrainLabel,Data.Test,Data.TestLabel,ParallelFlag,CSP_Config,1);
    [Acc(iFold,:), Kappa(iFold,:), tp{iFold}] = FBNN_SVM_Excute(PatFeature,accMat,SelFlag,1);
    t = toc;
    set(handles.TextState, 'String', ['Completed: ' num2str(iFold) '/' num2str(nFold) ' Fold']);
    disp(['Completed: ' num2str(iFold) '/' num2str(nFold) ' Fold']);
end

Acc_CV = mean(Acc);
Kappa_CV = mean(Kappa);

set(handles.TextResLinearSVM, 'String', [num2str(Acc_CV(1)*100,'%.2f'), '%']);
% set(handles.TextResQuadSVM, 'String', [num2str(Acc_CV(2)*100), '%']);
% set(handles.TextResCubeSVM, 'String', [num2str(Acc_CV(3)*100), '%']);
% set(handles.TextResBiquadSVM, 'String', [num2str(Acc_CV(4)*100), '%']);

% --- Executes on button press in ButtonBrowser.
function ButtonBrowser_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonBrowser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global TrainSetting

DataPath = uigetdir('*.*','Please select the folder');
DataPath = [DataPath, '\'];
set(handles.EditTrainDataPath, 'String', DataPath);
TrainSetting.TrainingDataFold = DataPath;

% disp(['数据加载完成，用时' num2str(toc) '秒']);

% --- Executes on button press in ButtonLoadData.
function ButtonLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global TrainSetting

nData = numel(dir(TrainSetting.TrainingDataFold)) - 2;
ImageData = [];
for iData = 1:nData
    load([TrainSetting.TrainingDataFold, num2str(iData), '.mat'], 'MotorImagery');
    TrainSetting.ImageData(iData,:,:) = MotorImagery.RecordData(:, MotorImagery.Sample.Stage2+1:MotorImagery.Sample.Stage2+1000);
    TrainSetting.ImageLabel(iData,:) = MotorImagery.ClassOrder;
end
% TrainSetting.ImageData = TrainSetting.ImageData(:,1:6,:);
if sum(mean(TrainSetting.ImageLabel)-TrainSetting.ImageLabel(1,:)) == 0
    TrainSetting.ImageLabel = MotorImagery.ClassOrder';
    set(handles.TextState, 'String', 'Data loaded successfully!');
%     disp(['数据标签验证成功！']);
else
    set(handles.TextState, 'String', 'Mislabeled data!');
    error('Mislabeled data!');
end


% --- Executes during object creation, after setting all properties.
function EditTrainParaFilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditTrainParaFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
