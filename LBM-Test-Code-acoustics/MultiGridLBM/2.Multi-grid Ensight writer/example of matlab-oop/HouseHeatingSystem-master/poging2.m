%%
% room<house
% 记录房间的温度
clc
%%
profile1 = [ 0,  6,  9, 17, 23, 24;
            16, 21, 16, 18, 16, 16];

          
%% Init
%实例化对象
woonkamer   = room('woonkamer', 20, 20);
slaapmaker1 = room('slaapkamer 1', 20, profile1);
slaapmaker2 = room('slaapkamer 2', 20, 20);
slaapmaker3 = room('slaapkamer 3', 20, 20);
serverkamer = room('server kamer', 10, 10);

% profile1 = [ 0,  6,  6,  9,  9, 17, 17, 23, 23, 24;
%             16, 16, 21, 21, 16, 16, 18, 18, 16, 16];

          
          
huisje = house([woonkamer, slaapmaker1, slaapmaker2, ...
                slaapmaker3, serverkamer]);
              
huisje.list_rooms;


%%
huisje.plot_a_day;

%%

huisje.plot_a_day_programming;




