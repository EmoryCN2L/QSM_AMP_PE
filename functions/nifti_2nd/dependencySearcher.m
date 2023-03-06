function [nameList,dep] = dependencySearcher(dirList,n)
%search for dependency

if(ischar(dirList))
    dirList = cellstr(dirList);
end;

dep = createDep(dirList);
id = [];
for i = 1:numberOfMat(dep);
    if(isequal(dep.itemList(i).name,n)||isequal([dep.itemList(i).name,'.m'],n))
        id = i;
        break;
    end;
end;

if(isempty(id))
    error('Can not find %',n);
end;

lst = getDepList(dep,id);

nameList = listToNameList(dep,lst);

return;




function dep = createDep(dirList)
    dep.itemList = zeros(1,0);%= struct('name','','fullpath','');
    %create a list of mat files
    for i = 1:length(dirList)
        lst = dir(fullfile(dirList{i},'*.m'));
        for j = 1:length(lst)
            dep.itemList(end+1).name=lst(j).name(1:(end-2));%delete .m extension
            dep.itemList(end).fullpath = fullfile(dirList{i},lst(j).name);
        end;
    end;
    numOfMat = numberOfMat(dep);
    dep.matrix = sparse(numOfMat,numOfMat);%initialize to be none dependency
   
    if(numberOfMat(dep)==0)
     error('no mat file found');
    end;

    for i = 1:numberOfMat(dep)
        str = readfile(dep.itemList(i).fullpath);
        for j = 1:numberOfMat(dep);
            if(~isempty(strfind(str',dep.itemList(j).name)))
                dep.matrix(i,j) = 1;
            end;
        end;

    end;
return;

function num = numberOfMat(dep)
    num = numel(dep.itemList);
return;

function d = readfile(n)
fid = fopen(n,'r');
d = fread(fid,inf,'char=>char');
fclose(fid);
return;

function lst = getDepList(dep,id,lst)
if(exist('lst','var')==0)
    lst=[];
end;

for j= 1:numberOfMat(dep)
    if(dep.matrix(id,j)==1&&isempty(find(lst==j,1)))
        lst(end+1) = j;%add current item to the list
        lst=getDepList(dep,j,lst);
    end;
end;

return;

function namelist = listToNameList(dep,lst)
namelist = {dep.itemList(lst).fullpath};
return;
