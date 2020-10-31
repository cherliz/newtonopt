function [xbest,fbest,itrcnt,stat] = NewtonOpt(f,Df,Df2,x0,gradTol,xTol,itrBound,s)
% INPUT VALUES:
% f - A filename or function handle that evaluates the function to be minimized.
% Df - A filename or function handle that evaluates the gradient of f.
% Df2 - A filename or function handle that evaluates the Hessian of f.
% x0 - An initial iterate.
% gradTol - A termination criterion on det(Dfx(xbest))
% xTol - A termination criterion on det(x(k)?x(k?1))
%        The algorithm should stop when either Dfx(xbest)
%        is sufficiently small or successive iterates 
%        are sufficiently close.
% itrBound - A maximum number of iterations to perform.
% s - A line search parameter encoded as follows:
%     s = 1: Use ? = 1 for all iterations.
%     s = 2: Decide ? with a Lagrange polynomial with ?1 = 0,
%            ?2 = 1, and ?3 = 2.
%     s = 3: Decide ? with bisection over the interval 0 ? ? ? 1;
%            if ?f(x) and ?f(x+d) have the same sign, then set
%            ? = 1.
%     s = 4: Decide ? by solving g'(s) = 0, where g(s) is the
%            second order Taylor polynomial for f(x + sd).
%
% OUTPUT VALUES:
%  xbest - The best calculated solution.
%  fxbest - The function value at xbest.
%  itrcnt - The number of iterations performed, with x0 counting as
%           iteration 0.
%  stat - A status variable, encoded as follows:
%         s = 0: The algorithm succeeded and  ?f(xbest)
%                is less than gradtol.
%         s = 1: The algorithm failed, in which case xbest is the best
%                solution that was encountered. 

xbest = x0;
% inv(A)*b can be less accurate than A/b

d  = -1 * (Df2(xbest) / Df(xbest));
a = calcA(s, xbest, d, f, Df, Df2);
i = 1;
while norm(Df(xbest))>=gradTol && i<=itrBound && norm(a*d)>=xTol
    d  = -1 * Df2(xbest) / Df(xbest);
    a = calcA(s, xbest, d, f, Df, Df2);
    xbest = xbest + a*d;
    i = i+a;
end

if norm(Df(xbest))<gradTol
   stat = 0;
else
    stat = 1;
end

itrcnt = i;
fbest = f(xbest);