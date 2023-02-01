%该函数随机初始化agent在搜索空间中的位置。
function [X]=initialization(SearchAgents_no,dim,up,down) %up,down两个参数随着目标函数的不同而对应不同的值，
%up,down实际上对应着蝗虫个体变量每一维的上下界，共dim维；给定目标函数后，给它们一个统一的上下界，大多数对应的上下界只是一个标量

if size(up,1)==1  % size(up,1):计算up的第一个维度的长度
    if SearchAgents_no>dim
        up(dim+1:SearchAgents_no)=up(1);
        down(dim+1:SearchAgents_no)=down(1);
    elseif SearchAgents_no<dim 
        up(SearchAgents_no+1:dim)=[];
        down(SearchAgents_no+1:dim)=[];
    end
    X=rand(SearchAgents_no,dim).*(up-down)+down;%得到元素位于[down,up]的N行dim列随机数矩阵(PS: .*代表对应元素相乘)
                                   %rand(N,dim)是指生成N行dim列的0-1之间的随机数构成的矩阵
   
end

if size(up,1)>1 %如果up的第一个维度的长度>1
    for i=1:dim
        high=up(i);
        low=down(i);
        X(:,i)=rand(1,SearchAgents_no).*(high-low)+low;
    end
end