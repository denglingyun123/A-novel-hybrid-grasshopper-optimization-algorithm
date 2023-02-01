function [x,f_x]=pattern_search1(x,fobj,dim,delta)

I=eye(dim);
alpha=1e-20;
a=0.5;
b=1;
n=dim;
y=x;
while min(abs(delta))>alpha
    for j=1:n
        f1=fobj(y+delta(j)*I(j,:));
        if f1<fobj(y)
            y=y+delta(j)*I(j,:);
        elseif fobj(y-delta(j)*I(j,:))<fobj(y)
            y=y-delta(j)*I(j,:);
        end
    end
    
    f3=fobj(y);
    if f3<fobj(x)
        x1=x;
        x=y;
        y=x+b*(x-x1);
    else 
        delta=delta*a;
        y=x;        
    end
end
f_x=fobj(x);




            


