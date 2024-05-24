clear all
close all
clc

X1=-10:0.01:10;
X2=-10:0.01:10;
[x1,x2]=meshgrid(X1,X2);
F=(1*cos((1+1)*x1+1)+2*cos((2+1)*x1+2)+3*cos((3+1)*x1+3)+4*cos((4+1)*x1+4)+5*cos((5+1)*x1+5))*(1*cos((1+1)*x2+1)+2*cos((2+1)*x2+2)+3*cos((3+1)*x2+3)+4*cos((4+1)*x2+4)+5*cos((5+1)*x2+5));
realFMin=min(min(F))
mesh(x1,x2,F) 

figure
contourf(x1,x2,F)
hold on
%% Newton-Raphson
fprintf('Newton-Raphson Algorithm\n');
% x=rand(2,1);
x0=-10+100*rand(2,1);
%x=[-1;1];
x=x0
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
plot(x(1),x(2),'r.')
x_next=x-inv(hessianfunc(x))*gradfunc(x);
fprintf('k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2),func(x_next),norm(gradfunc(x_next)))
%fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
plot(x_next(1),x_next(2),'r*')
k=3;
%while(abs(func(x_next)-func(x))>epsilon)
while(norm(gradfunc(x_next))>epsilon)
    x=x_next;
    x_next=x-inv(hessianfunc(x))*gradfunc(x);
    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n',k,x_next(1),x_next(2),func(x_next),norm(gradfunc(x_next)))
    %fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    plot(x_next(1),x_next(2),'r*')
    k=k+1;
end
toc
title('Newton-Raphson Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Hestenes-Stiefel Algorithm
figure
contourf(x1,x2,F)
hold on

fprintf('Hestenes-Stiefel Algorithm\n');
% x=rand(2,1);
%x=[-1;1];
x=x0;
epsilon=10^(-4);

tic
fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
plot(x(1),x(2),'r.')
g=gradfunc(x);
d=-g;

alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end
[val,ind]=min(funcalpha);
alpha=alpha(ind);

x_next=x+alpha*d;
g_next=gradfunc(x_next);
beta=(g_next'*(g_next-g))/(d'*(g_next-g));
d_next=-g_next+beta*d;

%fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
fprintf('k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2),func(x_next),norm(gradfunc(x_next)))
plot(x_next(1),x_next(2),'r*')
k=3;
while(norm(gradfunc(x_next))>epsilon)
    x=x_next;
    g=g_next;
    d=d_next;
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    
    x_next=x+alpha*d;
    g_next=gradfunc(x_next);
    beta=(g_next'*(g_next-g))/(d'*(g_next-g));
    d_next=-g_next+beta*d;

    %fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n',k,x_next(1),x_next(2),func(x_next),norm(gradfunc(x_next)))
    plot(x_next(1),x_next(2),'r*')
    k=k+1;
end
toc
title('Hestenes-Stiefel Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);