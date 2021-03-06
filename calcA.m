function [a] = calcA(s, x, d, f, Df ,Df2)

switch(s)
    case 1
        a = 1;
    case 2
        % decide alpha with a Lagrange polynomial with
        a1 = 0;
        a2 = 1;
        a3 = 2;
        x1 = x + a1*d;
        x2 = x + a2*d;
        x3 = x + a3*d;
        
        a = ((a2^2-a3^2)*f(x1) + (a3^2-a1^2)*f(x2) + (a1^2-a2^2)*f(x3)) ...
            /(2*(a2-a3)*f(x1) + 2*(a3-a1)*f(x2) + 2*(a1-a2)*f(x3));
    case 3
        % Bisection
        if ((Df(x)'*d<0) == (Df(x+d)'*d<0))
            a = 1;
        else
            itr = 50;
            l=0;
            h=1;
            
            for i=1:itr
                mid = (l+h)/2;
                if ((Df(x+mid)'*d<0)==(Df(x+d+mid)'*d<0))
                    h = mid;
                else
                    l = mid;
                end
            end
            a = mid;
        end
    case 4
        % solving g'(s)=0 where g(s) is second order taylor
        % poly for f(x+sd)
        a = -1*(Df(x+s*d)'*d)/(d'*(Df2(x+s*d))*d);
end