function data = UnpackData(data)
    dataFirst = data(:,1);
    data = cat(2, data(:,1)*256^2, data(:,2)*256, data(:,3));
    data = sum(data,2);
    data(dataFirst>=127) = data(dataFirst>=127)-256^3;
end