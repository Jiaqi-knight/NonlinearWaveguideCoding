function output=Q_lean(theta)
output=[  cos(theta)   sin(theta)  0;
          -sin(theta)   cos(theta)   0 ;
            0            0          1];%经过stagger变换后，A-4:绕z轴逆时针为正   
end