%�ú��������ʼ��agent�������ռ��е�λ�á�
function [X]=initialization(SearchAgents_no,dim,up,down) %up,down������������Ŀ�꺯���Ĳ�ͬ����Ӧ��ͬ��ֵ��
%up,downʵ���϶�Ӧ�Żȳ�������ÿһά�����½磬��dimά������Ŀ�꺯���󣬸�����һ��ͳһ�����½磬�������Ӧ�����½�ֻ��һ������

if size(up,1)==1  % size(up,1):����up�ĵ�һ��ά�ȵĳ���
    if SearchAgents_no>dim
        up(dim+1:SearchAgents_no)=up(1);
        down(dim+1:SearchAgents_no)=down(1);
    elseif SearchAgents_no<dim 
        up(SearchAgents_no+1:dim)=[];
        down(SearchAgents_no+1:dim)=[];
    end
    X=rand(SearchAgents_no,dim).*(up-down)+down;%�õ�Ԫ��λ��[down,up]��N��dim�����������(PS: .*�����ӦԪ�����)
                                   %rand(N,dim)��ָ����N��dim�е�0-1֮�����������ɵľ���
   
end

if size(up,1)>1 %���up�ĵ�һ��ά�ȵĳ���>1
    for i=1:dim
        high=up(i);
        low=down(i);
        X(:,i)=rand(1,SearchAgents_no).*(high-low)+low;
    end
end