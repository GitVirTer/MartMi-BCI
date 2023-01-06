function [net, sumh2, suma2, PredScore] = trainNet2_ButterFeature(Data, FinalTestData, DataTrainLabel, DataTestLabel)

sumh1 = []; sumh2 = [];
suma1 = []; suma2 = [];

for nbigcycle = 1:1
% clearvars -except nbigcycle sumh1 sumh2 suma1 suma2 acc1 acc2 Data FinalTestData indices FinalTestLabel DataTrainLabel DataTestLabel
% n_test = (indices == nbigcycle); %获得test集元素在数据集中对应的单元编号
% n_train = ~n_test;

% for iClass = 1:4
%     ClassMat(:,iClass) = find(DataTrainLabel == iClass);
% end
% nRow = round((1/2)*size(ClassMat, 1));
% TrainClassVector = ClassMat(1:nRow, :);
% TrainClassVector = TrainClassVector(:);
% TestClassVector = ClassMat(nRow+1:end, :);
% TestClassVector = TestClassVector(:);

XTrain = Data; 
% XTrain = Data(:,:,:, 1:nTrain);
% XTrain = Data(:,:,:, TrainClassVector);

YTrain = categorical(DataTrainLabel);
% YTrain = categorical(DataTrainLabel(1:nTrain));
% YTrain = categorical(DataTrainLabel(TrainClassVector));

% XTest = Data(:, :, :, nTrain+1:end);
XTest = FinalTestData;

% YTest = categorical(DataTrainLabel(nTrain+1:end));
YTest = categorical(DataTestLabel);

% FinalTestLabel = categorical(DataTestLabel);


numClasses = 2;

layers = [ ...
    imageInputLayer([size(XTrain,1) size(XTrain,2) size(XTrain,3)])

    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];
maxEpochs = 1000;
miniBatchSize = 300;

options = trainingOptions('adam', ...
    'InitialLearnRate',1e-2, ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'Shuffle','every-epoch', ...
    'Verbose',0, ...
    'LearnRateSchedule','none', ...
    'LearnRateDropFactor',0.05, ...
    'ExecutionEnvironment','cpu', ...
    'LearnRateDropPeriod',20);

[net, trainInfo] = trainNetwork(XTrain,YTrain,layers,options);

% [YPred,~] = classify(net,XTest);
miniBatchSize = 512;
[YPred_FinalTest,PredScore] = classify(net,XTest, ...
        'MiniBatchSize',miniBatchSize,...
        'ExecutionEnvironment','cpu');

% acc1(nbigcycle) = sum(YPred == YTest)./numel(YTest);
acc2(nbigcycle) = sum(YPred_FinalTest == YTest)./numel(YTest);

disp(['acc2: ', num2str(acc2(nbigcycle)*100), '%']);

% data = zeros(numel(res_test),numClasses);
% for n = 1:floor(numel(YTest)/numel(res_test))
%     data = data + YScore((n-1)*numel(res_test)+1:n*numel(res_test),:);
% end
% [~,dataRes] = max(data');



% sumh1 = [sumh1; double(YPred)];
% suma1 = [suma1; double(YTest)];

sumh2 = [sumh2; double(YPred_FinalTest)];
suma2 = [suma2; double(YTest)];
        
end