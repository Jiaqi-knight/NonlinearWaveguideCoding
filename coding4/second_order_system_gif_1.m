function second_order_system_gif_1()
    %此种情况用来表示某个参数变化时，对相关图形的影响
    %此函数用来表示当二阶系统的共轭极点的虚部变化时对系统的
    %频率响应和阶跃响应的影响
    %后期可以通过Ulead GIF Animator软件将这4张gif合并在一起
    %问题：图像大小设置为某些值时，可能会出错，需要重新调整
    
    clc;clear;close all;
    %初始化数据
    b = 2;
    a = 0:0.5:20;
    [~,size_a] = size(a);
    num = b.^2 + a.^2;
    for i = 1:size_a
        den{i} = conv([1 b + a(i) * 1i],[1  b - a(i) * 1i]);
    end
    w = 0:0.01:30;
    k = 1;
    %有多幅图，用胞元数组来指定文件名，从而方便在循环中使用
    fieldnames = {'1.gif','2.gif','3.gif','4.gif'};

    %画图，并制作gif
    %由于每次画图都要擦掉上一次画的图，所以图形不能一直用hold on
     for i = 1:1:size_a
        %完成图像的绘制，为了保证效果，要保证图像大小以及
        %坐标轴的范围不变
        figure(1);
         set(gcf,'Position',[0,0,300,400], 'color','w'); 
         
         [hz,hp,ht] = zplane(num(i),den{i});
         hold on;
         x_data = [0 hp.XData 0];
         y_data = [0 hp.YData 0];
         plot(x_data,y_data,'--');
         ylim([-22,22]);
         xlim([-6,6]);
         title(['二阶系统的极点',char(10,13)',...
             'Created by Lijunjie']);
         set(gca,'XTick',[-6:2:6]);
         hold off
        
        %采集绘制频率响应的数据
        h = freqs(num(i), den{i},w);
        mag = abs(h);
        phase = angle(h);
        phasedeg = phase*180/pi;

        figure(2)
        picture_positon;
        plot(w,mag);
        grid on;
        xlabel 'Frequency (rad/s)', ylabel Magnitude
        ylim([0 5.5]);
        xlim([0 30]);
        title(['二阶系统的幅频特性',char(10,13)',...
             'Created by Lijunjie']);
        figure(3);
        picture_positon;
        plot(w,phasedeg);
        xlabel 'Frequency (rad/s)', ylabel 'Phase (degrees)';
        ylim([-200,0]);
        xlim([0 30]);
        title(['二阶系统的相频特性',char(10,13)',...
             'Created by Lijunjie']); 
         
         sys = tf(num(i),den{i});
         figure(4)
         [y_tmp,t_tmp] = step(sys,3.5);
         plot(t_tmp,y_tmp);
         picture_positon;
         title(['二阶系统的阶跃响应',char(10,13)',...
             'Created by Lijunjie']); 
         xlabel('Time(seconds)');
         ylabel('Amplitude');
         axis([0 3.5 0 2]);
         
        %制作pdf
        if  i == 1
            %采集到首帧，需要设置gif的样式，以及确定图像的大小
            for j = 1:4
                figure(j)
                frame = getframe(gcf); % 获取整个窗口内容的图像
                im=frame2im(frame);
                [I{j,k},map{j,k}]=rgb2ind(im,256);
                imwrite(I{j,k},map{j,k},fieldnames{j},'gif','Loopcount',Inf,'DelayTime',0.2);
            end
        else
            for j = 1:4
                figure(j)
                frame = getframe(gcf);% 获取整个窗口内容的图像
                im=frame2im(frame);
                [I{j,k},map{j,k}]=rgb2ind(im,256);
                %追加模式
                imwrite(I{j,k},map{j,k},fieldnames{j},'gif','WriteMode','append','DelayTime',0.1);  
            end
        end
        k = k + 1;
     end
     
     %将采集到的图像以相反的顺序写入
     for i = (k-1):-1:1
         for j = 1:4
             imwrite(I{j,i},map{j,i},fieldnames{j},'gif','WriteMode','append','DelayTime',0.1); 
         end
     end
 
    function picture_positon
        %设置图像的大小
       set(gcf,'Position',[0,0,600,400], 'color','w'); 
    end
 end