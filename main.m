clear all
clc
close all
SearchAgents_no=30; % Number of search agents
Function_name='F1'; % Name of the test function 
Max_iteration=1000; % Maximum number of iterations  

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

Max_test=30;
for i=1:Max_test
    disp(['The ',num2str(i),'-th experiment']);
    [Best_pos2(i,:),Best_score2(i),HGOA_curve(i,:)]=HGOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
end

figure
semilogy(mean(HGOA_curve),'color','[0.62745,0.12549,0.94118]','linewidth',2.0,'Marker','o','MarkerIndices',1:100:length(mean(HGOA_curve)))
title('Convergence curve of F_{1}')
xlabel('Iteration');
ylabel('Fitness');
axis tight
grid off
box on 
legend('HGOA')


disp('-------------------------------------------------')
display(['The best fitness value obtained by HGOA over 30 independent runs : ', num2str(min(Best_score2))]);
display(['The mean fitness value obtained by HGOA over 30 independent runs : ', num2str(mean(Best_score2))]);
display(['The worst fitness value obtained by HGOA over 30 independent runs  : ', num2str(max(Best_score2))]);
display(['The standard deviation obtained by HGOA over 30 independent runs : ', num2str(std(Best_score2))]);

