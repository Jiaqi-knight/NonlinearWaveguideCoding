clear all
close all

%========================
% FSELIB
%
% Code: elm_line3_dr
%
% Element discretization
%=======================
                                                                                
%-----------
% input data
%-----------
                                                                                
x1 =-0.1;
x2 = 1.4;
n  = 21;        % must be odd
ratio = 10.0;
                                                                                
%-----------
% discretize
%-----------

xe = elm_line3 (x1,x2,n,ratio);

%-----
% plot
%-----

figure(1)
ye = zeros(n+1,1);
plot(xe, ye,'-ok');
set(gca,'fontsize',12)
axis equal
axis([x1-0.10,x2+0.10,-0.1,0.1])
axis off

%-----
% done
%-----
