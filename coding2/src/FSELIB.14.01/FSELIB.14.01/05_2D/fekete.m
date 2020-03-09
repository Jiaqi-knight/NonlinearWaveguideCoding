function [nodes,xi,et,xih,eth] = fekete(m)

%============================
% evaluates the Fekete points
% for m=6
%===========================

Ic=0;

if(m==6)

Ic=Ic+1;
xi(Ic) = 0.0;               et(Ic) = 0.0;          orbit(Ic) = 3;
Ic=Ic+1;
xi(Ic) = 0.0848854223;      et(Ic) = 0.0;          orbit(Ic) = 3;
Ic=Ic+1;
xi(Ic) = 0.2655651402;      et(Ic) = 0.0;          orbit(Ic) = 3;
Ic=Ic+1;
xi(Ic) = 0.5;               et(Ic) = 0.0;          orbit(Ic) = 3;
Ic=Ic+1;
xi(Ic) = 1.0-0.2655651402;  et(Ic) = 0.0;          orbit(Ic) = 3;
Ic=Ic+1;
xi(Ic) = 1.0-0.0848854223;  et(Ic) = 0.0;          orbit(Ic) = 3;
Ic=Ic+1;
xi(Ic) = 0.1063354684;      et(Ic) = 0.1063354684; orbit(Ic) = 3;
Ic=Ic+1;
xi(Ic) = 0.1171809171;      et(Ic) = 0.3162697959; orbit(Ic) = 6;
Ic=Ic+1;
xi(Ic) = 1.0/3.0;           et(Ic) = 1.0/3.0;      orbit(Ic) = 1;

end

%------------------------------------
% nodes over the equilateral triangle
% by orbiting
%------------------------------------

nd=Ic;

Jc = 0;
for i=1:nd

  xihh = xi(i)+0.5*et(i);
  ethh = 0.5*sqrt(3)*et(i);
  angle = 2.0*pi/orbit(i);
  xrot = 0.5;
  yrot = 0.5/sqrt(3);

  for j=1:orbit(i)
   cs = cos(j*angle);
   sn = sin(j*angle);
   Jc = Jc+1;
   xih(Jc) =  (xihh-xrot)*cs-(ethh-yrot)*sn + xrot;
   eth(Jc) =  (xihh-xrot)*sn+(ethh-yrot)*cs + yrot;
  end

end

nodes = Jc;
for i=1:nodes
 xi(i) = xih(i) - eth(i)/sqrt(3);
 et(i) = 2.0/sqrt(3)*eth(i);
end

%----
% done
%----

return
