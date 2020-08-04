
%% written by saeid zeynali- saedznl1374@gmail.com - msc student of power systems engineering at university of Tabriz - Iran

% this is a matlab function that calculates the lagrange multipliers of
% the maximum entropy problem (example provided)
% the enropy function is considered to be in form of p(x)=(lambda0+lambda1*x+lambda2*x.^2+lamda3*x^.3+lamda4*x.^4+lamda5*x^5+lamda6*x^6);
% which is ploynominal so the probibility density function would be in form of phi(x)=exp(p(x))
% with known lamda values obtained from this function

%% using instructions 
%xmin is minimum value of your probability density ,and xmax is the maximum
%of it. varargin (variable input argument) is the moments of your variables
%wich you want to find the problity density of it ,varargin must be at least
%two and maximum six. i'd like to give an exaxmple;

%% example: we have the first four non-central moments of a density function as m1,m2,m3,m4
%% and we want remake the function by maximum antripy method.

%fist write a function handle like --> 
%phi=@(x,la) exp(-la(1)-la(2)*x-la(3)*x.^2-la(4)*x.^3-la(5)*x.^4);
%[F,flag]=Maximum_Entropy(xmin,xmax,varargin) %varargin are the moments
%f=@(x) phi(x,F);
%f is the probibility density function ;

%% plot
% x=xmin:0.001:xmax;
% y=f(x);
% plot(x,y);
%% cautions 
% if you have 6 moments you should write phi=@(x,la) exp(-la(1)-la(2)*x-la(3)*x.^2-la(4)*x.^3-la(5)*x.^4-la(6)*x.^5-la(7)*x.^6);
% if problem is not converging increase nfe unde each case inside the
% function
% if its too slow decrese nfe under each case

function [F,flag,phi]=Maximum_Entropy(xmin,xmax,varargin) %varargin are the moments
      
         
         n=2500; %integration steps 
         
         sv=nargin-2;
        
         switch sv
             
             case 1
                 warning('only one moments is in the input, its not enough!')
                 
             case 2
                 nfe=30 ;%maximum number of function evaluation in fsolve
                 m1=varargin{1};
                 m2=varargin{2};
                 
                 phi=@(x,la) exp(-la(1)-la(2)*x-la(3)*x.^2);
                 ff0=@(x,la) phi(x,la);
                 f0=@(la) Trapezoid2(ff0,la,xmin,xmax,n)-1;

                 ff1=@(x,la) x*phi(x,la);
                 f1=@(la) Trapezoid2(ff1,la,xmin,xmax,n)-m1;

                 ff2=@(x,la) x^2*phi(x,la);
                 f2=@(la) Trapezoid2(ff2,la,xmin,xmax,n)-m2;
                 
                 FL=@(la) main2(la,f0,f1,f2)';
                 la0=[0.8 -0.3 0.8];
                 
             case 3
                 nfe=50; %maximum number of function evaluation in fsolve
                 m1=varargin{1};
                 m2=varargin{2};
                 m3=varargin{3};
                 
                 phi=@(x,la) exp(-la(1)-la(2)*x-la(3)*x.^2-la(4)*x.^3);
                 ff0=@(x,la) phi(x,la);
                 f0=@(la) Trapezoid2(ff0,la,xmin,xmax,n)-1;

                 ff1=@(x,la) x*phi(x,la);
                 f1=@(la) Trapezoid2(ff1,la,xmin,xmax,n)-m1;

                 ff2=@(x,la) x^2*phi(x,la);
                 f2=@(la) Trapezoid2(ff2,la,xmin,xmax,n)-m2;
                 
                 ff3=@(x,la) x^3*phi(x,la);
                 f3=@(la) Trapezoid2(ff3,la,xmin,xmax,n)-m3;
                 
                 FL=@(la) main3(la,f0,f1,f2,f3)';
                 la0=[0 0 0 0];
                 
             case 4
                 nfe=300; %maximum number of function evaluation in fsolve
                 m1=varargin{1};
                 m2=varargin{2};
                 m3=varargin{3};
                 m4=varargin{4};
                 
                 phi=@(x,la) exp(-la(1)-la(2)*x-la(3)*x.^2-la(4)*x.^3-la(5)*x.^4);
                 ff0=@(x,la) phi(x,la);
                 f0=@(la) Trapezoid2(ff0,la,xmin,xmax,n)-1;

                 ff1=@(x,la) x*phi(x,la);
                 f1=@(la) Trapezoid2(ff1,la,xmin,xmax,n)-m1;

                 ff2=@(x,la) x^2*phi(x,la);
                 f2=@(la) Trapezoid2(ff2,la,xmin,xmax,n)-m2;
                 
                 ff3=@(x,la) x^3*phi(x,la);
                 f3=@(la) Trapezoid2(ff3,la,xmin,xmax,n)-m3;
                 
                 ff4=@(x,la) x^4*phi(x,la);
                 f4=@(la) Trapezoid2(ff4,la,xmin,xmax,n)-m4;
                 
                 FL=@(la) main4(la,f0,f1,f2,f3,f4)';
                 la0=[0 0 0 0 0];
                 
             case 5
                 nfe=300; %maximum number of function evaluation in fsolve
                 m1=varargin{1};
                 m2=varargin{2};
                 m3=varargin{3};
                 m4=varargin{4};
                 m5=varargin{5};
                 
                 phi=@(x,la) exp(-la(1)-la(2)*x-la(3)*x.^2-la(4)*x.^3-la(5)*x.^4-la(6)*x.^5);
                 ff0=@(x,la) phi(x,la);
                 f0=@(la) Trapezoid2(ff0,la,xmin,xmax,n)-1;

                 ff1=@(x,la) x*phi(x,la);
                 f1=@(la) Trapezoid2(ff1,la,xmin,xmax,n)-m1;

                 ff2=@(x,la) x^2*phi(x,la);
                 f2=@(la) Trapezoid2(ff2,la,xmin,xmax,n)-m2;
                 
                 ff3=@(x,la) x^3*phi(x,la);
                 f3=@(la) Trapezoid2(ff3,la,xmin,xmax,n)-m3;
                 
                 ff4=@(x,la) x^4*phi(x,la);
                 f4=@(la) Trapezoid2(ff4,la,xmin,xmax,n)-m4;
                 
                 ff5=@(x,la) x^5*phi(x,la);
                 f5=@(la) Trapezoid2(ff5,la,xmin,xmax,n)-m5;
                 
                 
                 FL=@(la) main5(la,f0,f1,f2,f3,f4,f5)';
                 la0=[0 0 0 0 0 0];
                 
             case 6
                 nfe=400; %maximum number of function evaluation in fsolve
                 m1=varargin{1};
                 m2=varargin{2};
                 m3=varargin{3};
                 m4=varargin{4};
                 m5=varargin{5};
                 m6=varargin{6};
                 
                 phi=@(x,la) exp(-la(1)-la(2)*x-la(3)*x.^2-la(4)*x.^3-la(5)*x.^4-la(6)*x.^5-la(7)*x.^6);
                 ff0=@(x,la) phi(x,la);
                 f0=@(la) Trapezoid2(ff0,la,xmin,xmax,n)-1;

                 ff1=@(x,la) x*phi(x,la);
                 f1=@(la) Trapezoid2(ff1,la,xmin,xmax,n)-m1;

                 ff2=@(x,la) x^2*phi(x,la);
                 f2=@(la) Trapezoid2(ff2,la,xmin,xmax,n)-m2;
                 
                 ff3=@(x,la) x^3*phi(x,la);
                 f3=@(la) Trapezoid2(ff3,la,xmin,xmax,n)-m3;
                 
                 ff4=@(x,la) x^4*phi(x,la);
                 f4=@(la) Trapezoid2(ff4,la,xmin,xmax,n)-m4;
                 
                 ff5=@(x,la) x^5*phi(x,la);
                 f5=@(la) Trapezoid2(ff5,la,xmin,xmax,n)-m5;
                 
                 ff6=@(x,la) x^6*phi(x,la);
                 f6=@(la) Trapezoid2(ff6,la,xmin,xmax,n)-m6;
                 
                 
                 FL=@(la) main6(la,f0,f1,f2,f3,f4,f5,f6)';
                 la0=[0 0 0 0 0 0 0];
                 
             otherwise 
                 warning('too many input arguments! dont enter more than 6 moments.first input is xmin second is xmax and the rest are moments');
         end 
         options = optimoptions('fsolve','Display','off','MaxFunctionEvaluations',nfe,'FunctionTolerance',10e-6,'MaxIterations',80,'Algorithm','levenberg-marquardt');
         [L,fval,flag]=fsolve(FL,la0,options);
         
         display(['function value is:',num2str(sum(fval))]);
         
         if sum(fval)>0.06
             warning('Solver stopped prematurely, increase nfe!');
         end
         
         F=zeros(1,7);
         F(1:numel(L))=L;
         
end

function F=main2(la,f0,f1,f2)
F(1)=f0(la);
F(2)=f1(la);
F(3)=f2(la);
end

function F=main3(la,f0,f1,f2,f3)
F(1)=f0(la);
F(2)=f1(la);
F(3)=f2(la);
F(4)=f3(la);
end

function F=main4(la,f0,f1,f2,f3,f4)
F(1)=f0(la);
F(2)=f1(la);
F(3)=f2(la);
F(4)=f3(la);
F(5)=f4(la);
end

function F=main5(la,f0,f1,f2,f3,f4,f5)
F(1)=f0(la);
F(2)=f1(la);
F(3)=f2(la);
F(4)=f3(la);
F(5)=f4(la);
F(6)=f5(la);
end

function F=main6(la,f0,f1,f2,f3,f4,f5,f6)
F(1)=f0(la);
F(2)=f1(la);
F(3)=f2(la);
F(4)=f3(la);
F(5)=f4(la);
F(6)=f5(la);
F(7)=f6(la);
end

%% function to calculate integration
function s=Trapezoid2(f,la,a,b,n) %la is matrix; 
h=(b-a)/n;
s=f(a,la);

for i=a+h:h:b-h
    s=s+2*f(i,la);
end

s=(h/2)*(s+f(b,la));
end


%% flag
% 1
% Equation solved. First-order optimality is small.
% 2
% Equation solved. Change in x smaller than the specified tolerance.
% 3
% Equation solved. Change in residual smaller than the specified tolerance.
% 4
% Equation solved. Magnitude of search direction smaller than specified tolerance.
% 0
% Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations.
% -1
% Output function or plot function stopped the algorithm.
% -2
% Equation not solved. The exit message can have more information.
% -3
% Equation not solved. Trust region radius became too small (trust-region-dogleg algorithm).

        
        
        