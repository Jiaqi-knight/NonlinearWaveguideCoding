function h1=rotateIGV(h,azel,alpha,origin) % 旋转IGV并按照从前缘到尾缘的顺序输出，基准为0度


% find unit vector for axis of rotation
    u = azel(:)/norm(azel);


alph = (alpha+198.2)*pi/180;  %模型的角度为-18.2度
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = u(1);y = u(2);z = u(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';
   
   [m,n]=size(h); %newxyz = [x(:)-origin(1), y(:)-origin(2), z(:)-origin(3)];
   
   for k=1:n
    M=h{k};
    [m1,n1]=size(M);
    B=repmat(origin,m1,1);
    M1 = (M-B)*rot+B;
    M2=(fliplr(M1'))';%由于旋转接近180度，前尾缘倒置
    h1{k}=M2;   
   end
  
end
