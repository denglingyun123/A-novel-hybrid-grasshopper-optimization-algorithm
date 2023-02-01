%_________________________________________________________________________%
%  Grasshopper Optimization Algorithm (GOA) source codes demo V1.0        %
%                                                                         %
%  Developed in MATLAB R2016a                                             %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                                                                         %
%  Main paper: S. Saremi, S. Mirjalili, A. Lewis                          %
%              Grasshopper Optimisation Algorithm: Theory and Application %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2017.01.004   %
%                                                                         %
%_________________________________________________________________________%

% The Grasshopper Optimization Algorithm
function [TargetPosition,TargetFitness,Convergence_curve]=HGOA(N, Max_iter, lb,ub, dim, fobj)

disp('HGOA is now estimating the global optimum for your problem....')

flag=0;
if size(ub,1)==1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
end

if (rem(dim,2)~=0) % this algorithm should be run with a even number of variables. This line is to handle odd number of variables
    dim = dim+1;
    ub = [ub; 100];
    lb = [lb; -100];
    flag=1;
end

%Initialize the population of grasshoppers
GrassHopperPositions=initialization(N,dim,ub,lb);
GrassHopperFitness = zeros(1,N);

Convergence_curve=[];

cMax=1;
cMin=0.00004;
sensory_modality=0.01;
power_exponent=0.1;
%Calculate the fitness of initial grasshoppers

for i=1:size(GrassHopperPositions,1)
    if flag == 1
        GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,1:end-1));
    else
        GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,:));
    end
end

[sorted_fitness,sorted_indexes]=sort(GrassHopperFitness);

% Find the best grasshopper (target) in the first population 
for newindex=1:N
    Sorted_grasshopper(newindex,:)=GrassHopperPositions(sorted_indexes(newindex),:);
end

TargetPosition=Sorted_grasshopper(1,:);
TargetFitness=sorted_fitness(1);

Convergence_curve(1) = TargetFitness;

% Main loop
p=10^4;
l=2; % Start from the second iteration since the first iteration was dedicated to calculating the fitness of antlions
while l<=Max_iter
    tr=l/Max_iter;
    c=-p*cMin+(cMax+p*cMin)*exp(((tr).^1.5)*log((p+1)*cMin/(cMax+p*cMin)));
    if l<=Max_iter/2
        w = sin(pi*l/(Max_iter) + pi) + 1;
    else
        w = 0.15*cos(0.15*l*pi)+0.15;
    end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for i=1:size(GrassHopperPositions,1)
        temp= GrassHopperPositions';
       % for k=1:2:dim  
            S_i=zeros(dim,1);
            for j=1:N
                if i~=j
                    Dist=distance1(temp(:,j), temp(:,i)); % Calculate the distance between two grasshoppers
                    
                    r_ij_vec=(temp(:,j)-temp(:,i))/(Dist+eps); 
                    xj_xi=2+rem(Dist,2); 
                    
                    s_ij=((ub - lb)*c/2)*S_func(xj_xi).*r_ij_vec; % The first part inside the big bracket in Eq. (2.7)
                    S_i=S_i+s_ij;
                end
            end
            S_i_total = S_i;
            
      %  end
        
        X_new = c * S_i_total'+ (TargetPosition);    
        if rand<0.5
            f_new=fobj(X_new);
            if f_new<fobj(GrassHopperPositions(i,:))
                GrassHopperPositions_temp(i,:)=X_new';
            else
                GrassHopperPositions_temp(i,:)=GrassHopperPositions(i,:);
            end
        else
            Fnew=fobj(GrassHopperPositions(i,:));
            FP=(sensory_modality*(Fnew^power_exponent));   
            FP = real(FP);
            if rand>0.8
                Z1=w*GrassHopperPositions(i,:)+(rand*rand*TargetPosition-GrassHopperPositions(i,:))*FP;
                f_new1=fobj(Z1);
                if f_new1<fobj(GrassHopperPositions(i,:))
                    GrassHopperPositions_temp(i,:)=Z1;
                else
                    GrassHopperPositions_temp(i,:)=GrassHopperPositions(i,:);
                end
            else
                epsilon=rand;
                JK=randperm(N);
                Z2=w*GrassHopperPositions(i,:)+(epsilon*epsilon*GrassHopperPositions(JK(1),:)-GrassHopperPositions(JK(2),:))*FP;
                f_new2=fobj(Z2);
                if f_new2<fobj(GrassHopperPositions(i,:))
                    GrassHopperPositions_temp(i,:)=Z2;
                else
                    GrassHopperPositions_temp(i,:)=GrassHopperPositions(i,:);
                end
            end
        end
        
      end
      
     %% Centroid opposition-based learning
     sum3=0;
     for n=1:dim
         for m=1:N
             sum3=sum3+GrassHopperPositions_temp(m,n);
         end
         M(n)=sum3/N;
         sum3=0;
     end
     for m=1:N
         for n=1:dim
             OP(m,n)=2*M(n)-GrassHopperPositions_temp(m,n);
         end
         f_OP=fobj(OP(m,:));
         if f_OP<fobj(GrassHopperPositions_temp(m,:))
             GrassHopperPositions_temp(m,:)=OP(m,:);
         end
     end

    % GrassHopperPositions
    GrassHopperPositions=GrassHopperPositions_temp;
    
    for i=1:size(GrassHopperPositions,1)
        % Relocate grasshoppers that go outside the search space 
        Tp=GrassHopperPositions(i,:)>ub';
        Tm=GrassHopperPositions(i,:)<lb';
        GrassHopperPositions(i,:)=(GrassHopperPositions(i,:).*(~(Tp+Tm)))+ub'.*Tp+lb'.*Tm;
        
        % Calculating the objective values for all grasshoppers
        if flag == 1
            GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,1:end-1));
        else
            GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,:));
        end

        
        % Update the target
        if GrassHopperFitness(1,i)<TargetFitness
            TargetPosition=GrassHopperPositions(i,:);
            TargetFitness=GrassHopperFitness(1,i);
        end
    end
   
    [~,sorted_indexes3]=sort(GrassHopperFitness);
   M=(GrassHopperPositions(sorted_indexes3(1),:)+GrassHopperPositions(sorted_indexes3(N),:))/2;
   
   %% the step length in self-adaptive pattern search
    sum2=0;
    for j=1:dim
        for i=1:N
            sum2=sum2+GrassHopperPositions(i,j)-M(j);
        end
        delta(j)=sum2/N;
        sum2=0;
    end
    
    %% self-adaptive pattern search
    if rem(l,60)==0
        [TargetPosition1,~]=pattern_search1(TargetPosition,fobj,dim,delta);
        for a = 1:dim
            if(TargetPosition1(a)>ub(a))
                TargetPosition1(a) =ub(a);
            end
            if(TargetPosition1(a)<lb(a))
                TargetPosition1(a) =lb(a);
            end
       end

       f=fobj(TargetPosition1);
       if f<TargetFitness
           TargetPosition=TargetPosition1;
           TargetFitness=f;
       end
    end 
    
    sensory_modality=sensory_modality_NEW(sensory_modality,500);
   
    Convergence_curve(l)=TargetFitness;
    l=l+1;


%     disp(['Now the count of HGOA is:',num2str(TargetFitness)]);
end

if (flag==1)
    TargetPosition = TargetPosition(1:dim-1);
end


function y=sensory_modality_NEW(x,Ngen)
y=x+(0.025/(x*Ngen));



