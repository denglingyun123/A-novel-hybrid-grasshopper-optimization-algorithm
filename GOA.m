% The Grasshopper Optimization Algorithm
function [TargetPosition,TargetFitness,Convergence_curve]=GOA(N, Max_iter, lb,ub, dim, fobj)

disp('GOA is now estimating the global optimum for your problem....')

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

%Initialization
GrassHopperPositions=initialization(N,dim,ub,lb);
GrassHopperFitness = zeros(1,N);

Convergence_curve=[];

cMax=1;
cMin=0.00004;

for i=1:size(GrassHopperPositions,1)
    if flag == 1
        GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,1:end-1));
    else
        GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,:));
    end
end


[sorted_fitness,sorted_indexes]=sort(GrassHopperFitness);

% Find the best grasshopper in the first population 
for newindex=1:N
    Sorted_grasshopper(newindex,:)=GrassHopperPositions(sorted_indexes(newindex),:);
end

TargetPosition=Sorted_grasshopper(1,:);
TargetFitness=sorted_fitness(1);

Convergence_curve(1) = TargetFitness;

% Main loop
l=2; % Start from the second iteration since the first iteration was dedicated to calculating the fitness of antlions
while l<=Max_iter
    
    c=cMax-l*((cMax-cMin)/Max_iter); 
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for i=1:size(GrassHopperPositions,1)
        temp= GrassHopperPositions';
       % for k=1:2:dim  
            S_i=zeros(dim,1);
            for j=1:N
                if i~=j
                    Dist=distance1(temp(:,j), temp(:,i)); % Calculate the distance between two grasshoppers
                    
                    r_ij_vec=(temp(:,j)-temp(:,i))/(Dist+eps); % xj-xi/dij in Eq. (2.7)
                    xj_xi=2+rem(Dist,2); % |xjd - xid| in Eq. (2.7) 
                    
                    s_ij=((ub - lb)*c/2)*S_func(xj_xi).*r_ij_vec; % The first part inside the big bracket in Eq. (2.7)
                    S_i=S_i+s_ij;
                end
            end
            S_i_total = S_i;
            
      %  end
        
        X_new = c * S_i_total'+ (TargetPosition); % Eq. (2.7) in the paper      
        GrassHopperPositions_temp(i,:)=X_new'; 
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
    
     Convergence_curve(l)=TargetFitness;
     l=l+1; 
end


if (flag==1)
    TargetPosition = TargetPosition(1:dim-1);
end


