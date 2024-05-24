function ypr = gradfunc(x)
    syms x1 x2 
    F=(1*cos((1+1)*x1+1)+2*cos((2+1)*x1+2)+3*cos((3+1)*x1+3)+4*cos((4+1)*x1+4)+5*cos((5+1)*x1+5))*(1*cos((1+1)*x2+1)+2*cos((2+1)*x2+2)+3*cos((3+1)*x2+3)+4*cos((4+1)*x2+4)+5*cos((5+1)*x2+5));
    g=gradient(F);
    ypr=double(subs(g, [x1 x2], [x(1) x(2)]));
end