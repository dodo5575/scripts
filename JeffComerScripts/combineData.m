clear all;

nameList = {'0' '4' '8' '12' '20' '24' '30' '36' '39' '40' '41' '48' '52' '56' '60'};
prefix0 = 'str_Pz';
prefix1 = 'tilt_Pz';
outPrefix = 'versus';
suffix = '.dat';

for j=1:length(nameList)
    fileName0 = sprintf('%s%s%s', prefix0, nameList{j}, suffix);
    data0 = dlmread(fileName0, ' ');
    
    fileName1 = sprintf('%s%s%s', prefix1, nameList{j}, suffix);
    data1 = dlmread(fileName1, ' ');
    
    l = data0(:,2);
    a = data1(:,2);
    
    outName = sprintf('%s%s%s', outPrefix, nameList{j}, suffix);
    dlmwrite(outName ,[l a], ' ');
end
