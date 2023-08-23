%求目标函数的最大值
%测试基本完成
%程序初始化
clc
clear
close all

%%第1步，算法参数初始化
f=@(t) t.*sin(t);               %目标函数（性能指标）
t=0:0.01:30;
figure(1)
plot(t,f(t))
N = 50;                         %初始种群个数
d = 1;                          %空间维数（粒子的维度）
ger = 100;                      %最大迭代次数
xlimit = [0, 30];               %设置位置参数限制
vlimit = [-1, 1];               %设置速度限制
w = 0.8;                        %惯性权重
c1 = 0.5;                       %可信度因子1
c2 = 0.5;                       %可信度因子2
precision=0.0001;               %允许的精度误差
for i = 1:d
    x = xlimit(i, 1) + (xlimit(i, 2) - xlimit(i, 1)) * rand(N, d);%初始种群的位置
end
v = rand(N, d);                  % 初始种群的速度
Pbest = x;                          % 每个个体的历史最佳位置
Pbest_fit = zeros(N, 1);               % 每个个体的历史最佳适应度
for i = 1:N
    Pbest_fit(i) = f(x(i));        %初始种群每个粒子的适应度
end
Gbest = zeros(1, d);                % 种群的历史最佳位置
Gbest_fit = -inf;                      % 种群历史最佳适应度
hold on
plot(Pbest, f(Pbest), 'ro');
title('初始状态图');
figure(2)

%%开始迭代
iter = 1;
max_record = zeros(ger, 1);          % 记录每次迭代的最大值
while iter <= ger

    Gbestold=Gbest;
    %第2步， 计算每个粒子当前适应度
    fitness = f(x) ;
     %第3步，更新局部最优粒子和全局最优粒子，以及对应的适应度
     for i = 1:N
        if Pbest_fit(i) < fitness(i)
            Pbest_fit(i) = fitness(i);     % 更新个体历史最佳适应度
            Pbest(i,:) = x(i,:);   % 更新个体历史最佳位置
        end
     end
    if Gbest_fit < max(Pbest_fit)
        [Gbest_fit, nmax] = max(Pbest_fit);   % 更新群体历史最佳适应度
        Gbest = Pbest(nmax, :);      % 更新群体历史最佳位置
    end
    %第4步，速度更新和位置更新
    v = v * w + c1 * rand * (Pbest - x) + c2 * rand * (repmat(Gbest, N, 1) - x);% 速度更新
    v(v > vlimit(2)) = vlimit(2);% 速度边界速度处理
    v(v < vlimit(1)) = vlimit(1);

    x = x + v; % 位置更新
    x(x > xlimit(2)) = xlimit(2);% 位置边界处理
    x(x < xlimit(1)) = xlimit(1);
    %最大值记录
    max_record(iter) = Gbest_fit;
    x0 = 0 : 0.01 : 30;
    plot(x0, f(x0), 'b-', x, f(x), 'ro');title('状态位置变化')
    pause(0.1)


    %第5步，这里的判断加上外循环，判断迭代是否结束。满足收敛精度，同样可以终止迭代
    if(iter>50)
        if(abs(max_record(iter)-max_record(iter-1))<precision)&&(abs(Gbest-Gbestold)<precision)
            break
        end
    end
    iter = iter+1;



end
%最大值的收敛过程
figure(3);
plot(max_record(1:iter));
x0 = 0 : 0.01 : 30;
%粒子的收敛过程
figure(4);
plot(x0, f(x0), 'b-', x, f(x), 'ro'); 
