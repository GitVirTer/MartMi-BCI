function [Data_Train, Data_Test, Data_Train_Label, Data_Test_Label] = GetFoldData(DataTrain, DataTrainLabel, n_train, n_test, nClass)
    Data_Train = [];
    Data_Train_Label = [];
    Data_Test = [];
    Data_Test_Label = [];
    
    for iClass = 1:nClass
            DataClass = DataTrain(DataTrainLabel==iClass,:,:);
            
            DataClass_train = DataClass(n_train,:,:);
            Data_Train = cat(1,Data_Train,DataClass_train);
            DataClass_test = DataClass(n_test,:,:);
            Data_Test = cat(1,Data_Test,DataClass_test);            
            
            DataClassLabel = DataTrainLabel(DataTrainLabel==iClass,:);
            
            DataClassLabel_train = DataClassLabel(n_train,:);
            Data_Train_Label = cat(1,Data_Train_Label,DataClassLabel_train);
            DataClassLabel_test = DataClassLabel(n_test,:);
            Data_Test_Label = cat(1,Data_Test_Label,DataClassLabel_test);
            
    end
    
end