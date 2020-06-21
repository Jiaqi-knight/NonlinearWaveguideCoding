function output=Q_stagger(theta)
output=[ cos(theta)    0  -sin(theta) ;
          0            1     0
       sin(theta)     0  cos(theta);];%A-2:绕y轴逆时针为正
end