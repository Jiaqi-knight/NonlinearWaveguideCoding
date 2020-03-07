 t=linspace(0,2*pi,50);
tubeplot([cos(t);sin(t);0.2*(t-pi).^2],0.3);
 daspect([1,1,1]); camlight;