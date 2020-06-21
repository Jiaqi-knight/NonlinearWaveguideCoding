function output=Q_sweep(theta)
output=[  1    0             0 ;
          0  cos(theta)  sin(theta);
          0  -sin(theta)  cos(theta);];%经过stagger和lean变换后，A-6:绕x轴逆时针为正   
end