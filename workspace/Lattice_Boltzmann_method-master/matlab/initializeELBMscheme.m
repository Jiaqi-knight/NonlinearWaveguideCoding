function[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(D,Q,nongaussianflag)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 22nd, 2014
%    Last update: May 26th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

if nargin==2
    nongaussianflag = 0;
end

T0 = 0;
cs = sqrt(T0);
cssq = T0;
scheme = zeros(Q,D+2);

switch D
    case 1
        switch Q
            case 4
                m = 1;
                n = 4;
                Wm = ((m^2-5*n^2+sqrt(m^4-10*n^2*m^2+n^4))/(12*(m^2-n^2)));
                Wn = ((5*m^2-n^2-sqrt(m^4-10*n^2*m^2+n^4))/(12*(m^2-n^2)));
                T0 = (m^2+n^2+sqrt(m^4-10*n^2*m^2+n^4))/6;
                cs = sqrt(T0);
                cssq = T0;
                scheme = [m Wm 1;...
                          -m Wm 1;...
                          n Wn 4;...
                          -n Wn 4];
            case 5
                if nongaussianflag
                    W0 = 0.074464207985033;
                    W1 = 0.418585412256314;
                    W3 = 0.044182483751169;
                    T0 = 1.632455532033676;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme = [0 W0 0;...
                              1 W1 1;...
                             -1 W1 1;...
                              3 W3 3;...
                             -3 W3 3];
                else
                    m = 1;
                    n = 3;
                    W0 = (-3*m^4-3*n^4+54*m^2*n^2+(m^2+n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(75*m^2*n^2);
                    Wm = ((9*m^4-6*n^4-27*m^2*n^2-(3*m^2-2*n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*m^2*(m^2-n^2)));
                    Wn = ((9*n^4-6*m^4-27*m^2*n^2-(3*n^2-2*m^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*n^2*(n^2-m^2)));
                    T0 = (3*m^2+3*n^2-sqrt(9*m^4-42*n^2*m^2+9*n^4))/30;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme = [0 W0 0;...
                              m Wm m;...
                             -m Wm m;...
                              n Wn n;...
                             -n Wn n];
                end
            case 7
                m = 2;
                n = 3;
                % general formulation
                msq = m^2;
                nsq = n^2;
                mfp = m^4;
                nfp = n^4;
                r = n^2+1;
                a = m^6;
                b = r*mfp;
                c = (11*nfp-104*nsq+11)*msq;
                d = r*(10*nfp-43*nsq+10);
                h = msq + nsq + 1;
                e = 4*(945*(r*msq+nsq)-225*(h^2))^3;
                f = 455625*(10*a-33*r*mfp+(-33*nfp+312*nsq-33)*msq+r*(10*nfp-43*nsq+10))^2;
                g = e + f;
                i = r*msq + nsq; 
                l = msq*nsq;
                o = msq+nsq;
                q = msq+1;
                alpha = -250*a +825*b +75*c -25*d +(1/27)*sqrt(g);
                beta = -50*a +165*b +15*c -5*d +(1/135)*sqrt(g);
                gamma = 5*mfp-11*r*msq+5*nfp-11*nsq+5;
                %
                T0 = (10*h-(2^(2/3))*(alpha^(1/3))-2*(2^(1/3))*(5^(2/3))*gamma*(beta^(-1/3)))/210;
                %
                T0sq = T0^2;
                %
                T0tp = T0^3;
                %
                s = -15*T0tp+3*h*T0sq-i*T0+l;
                %
                W0 = s/l;
                %
                %
                t = (T0-1)*msq-3*T0sq+T0;
                %
                den = 2*(msq-1)*t;
                %
                aleph = l+3*T0sq-o*T0;
                %
                beth = msq+15*T0sq-3*q*T0;
                %
                gimel = -15*T0tp+3*h*T0sq-i*T0+l;
                %
                daleth = (T0-1)*nsq-3*T0sq+T0;
                %
                he = nsq+15*T0sq-3*r*T0;
                %
                tet = daleth*msq+T0*he;
                %
                firstaddendW1 = (T0*aleph*beth*gimel)/(nsq*(nsq-1)*tet);
                secondaddendW1 = ((msq-3*T0)*T0*s)/nsq;
                W1 = (firstaddendW1-secondaddendW1)/den;
                %
                firstaddendWm = (T0*(3*T0-1)*s)/l;
                chet = msq+15*T0sq-3*q*T0;
                waw = T0*t*chet*s;
                zajin = (m-n)*nsq*(m+n)*tet;
                secondaddendWm = waw/zajin;
                %
                Wm = -(firstaddendWm+secondaddendWm)/den;
                %
                Wn = (T0*chet*s)/(2*(m-n)*nsq*(m+n)*(nsq-1)*tet);
                %
                cs = sqrt(T0);
                cssq = T0;
                scheme = [0 W0 0;...
                          1 W1 1;...
                          -1 W1 1;...
                          m Wm m;...
                          -m Wm m;...
                          n Wn n;...
                          -n Wn n];
            case 9
                if nongaussianflag
                    T0 = 2.175382386573002;
                    W0 = 0.167240236271550;
                    W1 = 0.303154150433653;
                    W2 = 0.053302940534918;
                    W3 = 0.057921530729887;
                    W5 = 0.002001260165768;
                    scheme = [0 W0 0;...
                              1 W1 1;...
                             -1 W1 1;...
                              2 W2 2;...
                             -2 W2 2;...
                              3 W3 3;...
                             -3 W3 3;...
                              5 W5 5;...
                             -5 W5 5];
                else
                    m = 5;
                    p = m;
                    psq = p^2;
                    T0 =  0.756080852595360;
                    %
                    denW0 = 36*psq;
                    coeffT0W0 = 3*(7*T0*(5*(T0-2)*T0+7)-12);
                    coeffpsqW0 = T0*(3*(14-5*T0)*T0-49)+36;
                    %
                    W0 = (coeffpsqW0*psq+coeffT0W0*T0)/denW0;
                    %
                    denW1 = 16*(psq-1);
                    coeffT0W1 = 5*(13-7*T0)*T0-36;
                    coeffpsqW1 = (5*T0-13)*T0+12;
                    %
                    W1 = T0*(coeffpsqW1*psq+coeffT0W1*T0)/denW1;
                    %
                    denW2 = 40*(psq-4);
                    coeffT0W2 = 5*(7*T0-10)*T0+9;
                    coeffpsqW2 = -5*(T0-2)*T0-3;
                    %
                    W2 = T0*(coeffpsqW2*psq+coeffT0W2*T0)/denW2;
                    %
                    denW3 = 720*(psq-9);
                    coeffT0W3 = -3*(5*(7*T0-5)*T0+4);
                    coeffpsqW3 = 15*(T0-1)*T0+4;
                    %
                    W3 = T0*(coeffpsqW3*psq+coeffT0W3*T0)/denW3;
                    %
                    denWp = 2*psq*(psq*(psq-7)^2-36);
                    numWp = 7*T0*(5*(T0-2)*T0+7)-12;
                    %
                    Wp = 3*T0*numWp/denWp;
                    %
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme = [0 W0 0;...
                              1 W1 1;...
                             -1 W1 1;...
                              2 W2 2;...
                             -2 W2 2;...
                              3 W3 3;...
                             -3 W3 3;...
                              p Wp p;...
                             -p Wp p];
                end
        end       
    case 2
        switch Q
            case 25
                if nongaussianflag
                    W0 = 0.074464207985033;
                    W1 = 0.418585412256314;
                    W3 = 0.044182483751169;
                    T0 = 1.632455532033676;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              1 W1 1;...
                             -1 W1 1;...
                              3 W3 3;...
                             -3 W3 3];
                    [scheme,QD25] = constructentropicDnscheme(2,scheme1D);
                else
                    m = 1;
                    n = 3;
                    W0 = (-3*m^4-3*n^4+54*m^2*n^2+(m^2+n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(75*m^2*n^2);
                    Wm = ((9*m^4-6*n^4-27*m^2*n^2-(3*m^2-2*n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*m^2*(m^2-n^2)));
                    Wn = ((9*n^4-6*m^4-27*m^2*n^2-(3*n^2-2*m^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*n^2*(n^2-m^2)));
                    T0 = (3*m^2+3*n^2-sqrt(9*m^4-42*n^2*m^2+9*n^4))/30;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              m Wm m;...
                             -m Wm m;...
                              n Wn n;...
                             -n Wn n];
                    [scheme,QD25] = constructentropicDnscheme(2,scheme1D);
                end
                
        end
    case 3
        switch Q
            case 15
                if nongaussianflag
                    W0 = 0.074464207985033;
                    W1 = 0.418585412256314;
                    W3 = 0.044182483751169;
                    T0 = 1.632455532033676;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              1 W1 1;...
                             -1 W1 1;...
                              3 W3 3;...
                             -3 W3 3];
                    pruneflag = 1;
                    prune = [2;9;10;11;18;19;27];
                    [scheme,QD15] = constructentropicDnscheme(3,scheme1D,pruneflag,prune);
                else
                    m = 1;
                    n = 3;
                    W0 = (-3*m^4-3*n^4+54*m^2*n^2+(m^2+n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(75*m^2*n^2);
                    Wm = ((9*m^4-6*n^4-27*m^2*n^2-(3*m^2-2*n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*m^2*(m^2-n^2)));
                    Wn = ((9*n^4-6*m^4-27*m^2*n^2-(3*n^2-2*m^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*n^2*(n^2-m^2)));
                    T0 = (3*m^2+3*n^2-sqrt(9*m^4-42*n^2*m^2+9*n^4))/30;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              m Wm m;...
                             -m Wm m;...
                              n Wn n;...
                             -n Wn n];
                    pruneflag = 1;
                    prune = [2;9;10;11;18;19;27];
                    [scheme,QD15] = constructentropicDnscheme(3,scheme1D,pruneflag,prune);
                end
            case 19
                if nongaussianflag
                    W0 = 0.074464207985033;
                    W1 = 0.418585412256314;
                    W3 = 0.044182483751169;
                    T0 = 1.632455532033676;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              1 W1 1;...
                             -1 W1 1;...
                              3 W3 3;...
                             -3 W3 3];
                    pruneflag = 1;
                    prune = [3;9;10;11;18;19;27];
                    [scheme,QD15] = constructentropicDnscheme(3,scheme1D,pruneflag,prune);
                else
                    m = 1;
                    n = 3;
                    W0 = (-3*m^4-3*n^4+54*m^2*n^2+(m^2+n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(75*m^2*n^2);
                    Wm = ((9*m^4-6*n^4-27*m^2*n^2-(3*m^2-2*n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*m^2*(m^2-n^2)));
                    Wn = ((9*n^4-6*m^4-27*m^2*n^2-(3*n^2-2*m^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*n^2*(n^2-m^2)));
                    T0 = (3*m^2+3*n^2-sqrt(9*m^4-42*n^2*m^2+9*n^4))/30;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              m Wm m;...
                             -m Wm m;...
                              n Wn n;...
                             -n Wn n];
                    pruneflag = 1;
                    prune = [3;9;10;11;18;19;27];
                    [scheme,QD15] = constructentropicDnscheme(3,scheme1D,pruneflag,prune);
                end
            case 27
                if nongaussianflag
                    W0 = 0.074464207985033;
                    W1 = 0.418585412256314;
                    W3 = 0.044182483751169;
                    T0 = 1.632455532033676;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              1 W1 1;...
                             -1 W1 1;...
                              3 W3 3;...
                             -3 W3 3];
                    pruneflag = 1;
                    prune = [9;10;11;18;19;27];
                    [scheme,QD27] = constructentropicDnscheme(3,scheme1D,pruneflag,prune);
                else
                    m = 1;
                    n = 3;
                    W0 = (-3*m^4-3*n^4+54*m^2*n^2+(m^2+n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(75*m^2*n^2);
                    Wm = ((9*m^4-6*n^4-27*m^2*n^2-(3*m^2-2*n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*m^2*(m^2-n^2)));
                    Wn = ((9*n^4-6*m^4-27*m^2*n^2-(3*n^2-2*m^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*n^2*(n^2-m^2)));
                    T0 = (3*m^2+3*n^2-sqrt(9*m^4-42*n^2*m^2+9*n^4))/30;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              m Wm m;...
                             -m Wm m;...
                              n Wn n;...
                             -n Wn n];
                    pruneflag = 1;
                    prune = [9;10;11;18;19;27];
                    [scheme,QD27] = constructentropicDnscheme(3,scheme1D,pruneflag,prune);
                end
            case 41
                if nongaussianflag
                    W0 = 0.074464207985033;
                    W1 = 0.418585412256314;
                    W3 = 0.044182483751169;
                    T0 = 1.632455532033676;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              1 W1 1;...
                             -1 W1 1;...
                              3 W3 3;...
                             -3 W3 3];
                    pruneflag = 1;
                    prune = [10;11;18;19];
                    [scheme,QD41] = constructentropicDnscheme(3,scheme1D,pruneflag,prune);
                else
                    m = 1;
                    n = 3;
                    W0 = (-3*m^4-3*n^4+54*m^2*n^2+(m^2+n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(75*m^2*n^2);
                    Wm = ((9*m^4-6*n^4-27*m^2*n^2-(3*m^2-2*n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*m^2*(m^2-n^2)));
                    Wn = ((9*n^4-6*m^4-27*m^2*n^2-(3*n^2-2*m^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*n^2*(n^2-m^2)));
                    T0 = (3*m^2+3*n^2-sqrt(9*m^4-42*n^2*m^2+9*n^4))/30;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              m Wm m;...
                             -m Wm m;...
                              n Wn n;...
                             -n Wn n];
                    pruneflag = 1;
                    prune = [10;11;18;19];
                    [scheme,QD41] = constructentropicDnscheme(3,scheme1D,pruneflag,prune);
                    W00 = 2*(5045-1507*sqrt(10))/2025;
                    W100 = 37/(5*sqrt(10))-91/40;
                    W110 = (55-17*sqrt(10))/50;
                    W111 = (233*sqrt(10)-730)/1600;
                    W300 = (295-92*sqrt(10))/16200;
                    W333 = (130-41*sqrt(10))/129600;
                    weights = [W00*ones(1,1);W100*ones(6,1);W110*ones(12,1);W111*ones(8,1);W300*ones(6,1);W333*ones(8,1)];
                    scheme(:,end-1) = weights;
                end
            case 125
                if nongaussianflag
                    W0 = 0.074464207985033;
                    W1 = 0.418585412256314;
                    W3 = 0.044182483751169;
                    T0 = 1.632455532033676;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              1 W1 1;...
                             -1 W1 1;...
                              3 W3 3;...
                             -3 W3 3];
                    [scheme,QD125] = constructentropicDnscheme(3,scheme1D);
                else
                    m = 1;
                    n = 3;
                    W0 = (-3*m^4-3*n^4+54*m^2*n^2+(m^2+n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(75*m^2*n^2);
                    Wm = ((9*m^4-6*n^4-27*m^2*n^2-(3*m^2-2*n^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*m^2*(m^2-n^2)));
                    Wn = ((9*n^4-6*m^4-27*m^2*n^2-(3*n^2-2*m^2)*sqrt(9*m^4-42*n^2*m^2+9*n^4))/(300*n^2*(n^2-m^2)));
                    T0 = (3*m^2+3*n^2-sqrt(9*m^4-42*n^2*m^2+9*n^4))/30;
                    cs = sqrt(T0);
                    cssq = T0;
                    scheme1D = [0 W0 0;...
                              m Wm m;...
                             -m Wm m;...
                              n Wn n;...
                             -n Wn n];
                    [scheme,QD125] = constructentropicDnscheme(3,scheme1D);
                end
        end
end

invcs = 1/cs;
invcssq = 1/cssq;

return