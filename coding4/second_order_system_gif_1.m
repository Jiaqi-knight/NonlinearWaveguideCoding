function second_order_system_gif_1()
    %�������������ʾĳ�������仯ʱ�������ͼ�ε�Ӱ��
    %�˺���������ʾ������ϵͳ�Ĺ������鲿�仯ʱ��ϵͳ��
    %Ƶ����Ӧ�ͽ�Ծ��Ӧ��Ӱ��
    %���ڿ���ͨ��Ulead GIF Animator�������4��gif�ϲ���һ��
    %���⣺ͼ���С����ΪĳЩֵʱ�����ܻ������Ҫ���µ���
    
    clc;clear;close all;
    %��ʼ������
    b = 2;
    a = 0:0.5:20;
    [~,size_a] = size(a);
    num = b.^2 + a.^2;
    for i = 1:size_a
        den{i} = conv([1 b + a(i) * 1i],[1  b - a(i) * 1i]);
    end
    w = 0:0.01:30;
    k = 1;
    %�ж��ͼ���ð�Ԫ������ָ���ļ������Ӷ�������ѭ����ʹ��
    fieldnames = {'1.gif','2.gif','3.gif','4.gif'};

    %��ͼ��������gif
    %����ÿ�λ�ͼ��Ҫ������һ�λ���ͼ������ͼ�β���һֱ��hold on
     for i = 1:1:size_a
        %���ͼ��Ļ��ƣ�Ϊ�˱�֤Ч����Ҫ��֤ͼ���С�Լ�
        %������ķ�Χ����
        figure(1);
         set(gcf,'Position',[0,0,300,400], 'color','w'); 
         
         [hz,hp,ht] = zplane(num(i),den{i});
         hold on;
         x_data = [0 hp.XData 0];
         y_data = [0 hp.YData 0];
         plot(x_data,y_data,'--');
         ylim([-22,22]);
         xlim([-6,6]);
         title(['����ϵͳ�ļ���',char(10,13)',...
             'Created by Lijunjie']);
         set(gca,'XTick',[-6:2:6]);
         hold off
        
        %�ɼ�����Ƶ����Ӧ������
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
        title(['����ϵͳ�ķ�Ƶ����',char(10,13)',...
             'Created by Lijunjie']);
        figure(3);
        picture_positon;
        plot(w,phasedeg);
        xlabel 'Frequency (rad/s)', ylabel 'Phase (degrees)';
        ylim([-200,0]);
        xlim([0 30]);
        title(['����ϵͳ����Ƶ����',char(10,13)',...
             'Created by Lijunjie']); 
         
         sys = tf(num(i),den{i});
         figure(4)
         [y_tmp,t_tmp] = step(sys,3.5);
         plot(t_tmp,y_tmp);
         picture_positon;
         title(['����ϵͳ�Ľ�Ծ��Ӧ',char(10,13)',...
             'Created by Lijunjie']); 
         xlabel('Time(seconds)');
         ylabel('Amplitude');
         axis([0 3.5 0 2]);
         
        %����pdf
        if  i == 1
            %�ɼ�����֡����Ҫ����gif����ʽ���Լ�ȷ��ͼ��Ĵ�С
            for j = 1:4
                figure(j)
                frame = getframe(gcf); % ��ȡ�����������ݵ�ͼ��
                im=frame2im(frame);
                [I{j,k},map{j,k}]=rgb2ind(im,256);
                imwrite(I{j,k},map{j,k},fieldnames{j},'gif','Loopcount',Inf,'DelayTime',0.2);
            end
        else
            for j = 1:4
                figure(j)
                frame = getframe(gcf);% ��ȡ�����������ݵ�ͼ��
                im=frame2im(frame);
                [I{j,k},map{j,k}]=rgb2ind(im,256);
                %׷��ģʽ
                imwrite(I{j,k},map{j,k},fieldnames{j},'gif','WriteMode','append','DelayTime',0.1);  
            end
        end
        k = k + 1;
     end
     
     %���ɼ�����ͼ�����෴��˳��д��
     for i = (k-1):-1:1
         for j = 1:4
             imwrite(I{j,i},map{j,i},fieldnames{j},'gif','WriteMode','append','DelayTime',0.1); 
         end
     end
 
    function picture_positon
        %����ͼ��Ĵ�С
       set(gcf,'Position',[0,0,600,400], 'color','w'); 
    end
 end