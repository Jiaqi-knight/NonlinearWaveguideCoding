function H=plotCellData(h,color)
   [m,n]=size(h);
   for k=1:n
    M=h{k};
    H(k)=plot3(M(:,1),M(:,2),M(:,3),color);hold on
   end
end 
