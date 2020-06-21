function [ Arrow ] = arrow3d_relative(orgin,h1,r1,r2)
    Arrow(1)=arrow(orgin,orgin+[h1; 0; 0],r1,'BaseAngle',45,'Width',r2,'Color','g');
    Arrow(2)=arrow(orgin,orgin+[0; h1+40; 0],r1,'BaseAngle',45,'Width',r2,'Color','g');
    Arrow(3)=arrow(orgin,orgin+[0; 0; h1],r1,'BaseAngle',45,'Width',r2,'Color','g');
end
