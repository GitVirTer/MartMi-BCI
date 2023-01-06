function [nn, sumh2, suma2, PredScore, numB, Mdl1] = trainSVM_2Class(XTrain, XTest, YTrain, YTest, SelClass)
% 'Lambda',0.01, 'Alpha', 1   THR:0.3 %0.63
% 'Lambda',0.1, 'Alpha', 0.5   THR:0.2 %0.61
nn = {};sumh2 = 0; suma2 = 0;


if nargin > 4
    YTrain(YTrain==SelClass) = 0.2;
    YTrain(YTrain~=0.2) = -0.2;
else
    YTrain(YTrain==min(YTrain)) = -0.2;
    YTrain(YTrain==max(YTrain)) = 0.2;
end

YTrain(YTrain==-0.2) = 1;
YTrain(YTrain==0.2) = 2;
YTrain = categorical(YTrain);

Mdl1 = fitcsvm(XTrain, YTrain, 'KernelFunction', 'linear');

if ~isempty(XTest)
    Mdl2 = fitcsvm(XTrain, YTrain, 'KernelFunction', 'polynomial', 'PolynomialOrder', 2);
    Mdl3 = fitcsvm(XTrain, YTrain, 'KernelFunction', 'polynomial', 'PolynomialOrder', 3);
    Mdl4 = fitcsvm(XTrain, YTrain, 'KernelFunction', 'polynomial', 'PolynomialOrder', 4);

    [~,PredScore{1}] = predict(Mdl1, XTest);
    [~,PredScore{2}] = predict(Mdl2, XTest);
    [~,PredScore{3}] = predict(Mdl3, XTest);
    [~,PredScore{4}] = predict(Mdl4, XTest);
    % [~,PredScore{5}] = predict(Mdl5, XTest);
    % [~,PredScore{6}] = predict(Mdl6, XTest);

    numB = numel(PredScore);
else
    PredScore = [];
    numB = 0;
end

% % LambdaValue = [0.0001:0.0001:0.0009 0.001:0.001:0.009 0.01:0.01:0.09 0.1:0.05:0.8];
% [B,FitInfo] = lasso(XTrain, YTrain, 'Lambda', LambdaValue, 'Alpha', 0.07);
% for i = 1:size(B,2) data = B(:,i); Bnot0(i) = numel(data(data~=0)); end
% disp(['Not zero number: ' num2str(Bnot0)]);
% 
% % PredScore = zeros(size(XTest,1),2);
% for iB = 1:size(B,2)
%     LR = sigm(XTest*B(:,iB)+FitInfo.Intercept(iB));
%     PredScore{iB}(:,1) = 1-LR;
%     PredScore{iB}(:,2) = LR;
% end



% function acc = lassoCV(XTrain, YTrain, XTest, YTest)
%     YTrain(YTrain==min(YTrain)) = -0.2;
%     YTrain(YTrain==max(YTrain)) = 0.2;
% 
%     [B,FitInfo] = lasso(XTrain, YTrain, 'Lambda',[0.0001:0.0001:0.0009...
%                                                   0.001:0.001:0.009...
%                                                   0.01:0.01:0.09...
%                                                   0.1:0.05:0.3], 'Alpha', 0.5);
%     for i = 1:size(B,2) data = B(:,i); Bnot0(i) = numel(data(data~=0)); end
%     disp(['Not zero number: ' num2str(Bnot0)]);
% 
%     for iB = 1:size(B,2)
%         LR = sigm(XTest*B(:,iB)+FitInfo.Intercept(iB));
% %       LR = XTest*B+FitInfo.Intercept;
% 
%         Score(:,1) = 1-LR;
%         Score(:,2) = LR;
%         [~,res] = max(Score.');
%         acc(iB) = numel(find(res'==YTest))/numel(YTest);
%     end
%     
% end

end


