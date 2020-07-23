function [x,y,z,ds,dt]=pipesetup()

  a = 1; %the parameter 'a' from the equations, here a=2
%For the cylinder:
  t=linspace(0,2*pi,51); %defines the scope of 't'
  z=linspace(-2*a,2*a,30); %defines the scope of 'z', the bounds for the top and bottom of the cylinder
  t(end)=[];
  dt=2*pi/50;
  ds=4*a/30;

  [t,z]=meshgrid(t,z); %creates pairs of (t,z)
  x=a*cos(t); %the parametrized cylinder
  y=a*sin(t);
  z = 1*z;
%Plotting the cylinder
  hold on %prevents override of plot data, so both the cylinder and sphere can occupy the same plot
  mhndl2=mesh(x, y, z); %notice 'mhndl' is numbered, mhndl2=mesh(x,y,z) plots the cylinder
  set(mhndl2,...
  'EdgeColor',[1,0,0]) %[1,0,0] gives the color red
  view(50,20)
  axis equal
%view(horizontal_rotation,vertical_elevation) sets the angle of view,
%both parameters used are measured in degrees
end