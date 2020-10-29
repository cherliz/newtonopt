
type = 'confirmed'; % 'confirmed','deaths','recovered'

[dataMatrix] = readCoronaData(type);

[dataTable,timeVector,mergedData] = processCoronaData(dataMatrix);

usData = cell2mat(mergedData(176,2));

% plotCoronaData(timeVector,mergedData,{'Denmark','US','Germany','China','Italy'},type);

plot(timeVector, usData)
