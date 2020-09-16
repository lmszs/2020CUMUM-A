xls = xlsread('附件.xlsx','A2:B710');
tim = xls(:,1);
t = xls(:,2);
v = 7/6*0.01;
len = 435.5e-2;
% plot(tim,t);
KK = 273.15;

ll = 15e-2;    % 宽度
h = 0.5 / 100;
k = 0.1 / 100;
x = 0:h:(len);
y = 0:k:(ll);
%% t_c灵敏度分析
t_e = xlsread('z0.xls') + KK;
tc_real = 48;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
t_0(1) = 25 + KK;
for i = 2:length(ti)
    t_f = t_e(get_t2(ti(i),v,x));
    t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
end
plot(ti,t_0-KK);
hold on
t_e = xlsread('z0.xls') + KK;
tc_real = 43;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
t_0(1) = 25 + KK;
for i = 2:length(ti)
    t_f = t_e(get_t2(ti(i),v,x));
    t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
end
plot(ti,t_0-KK);
hold on
t_e = xlsread('z0.xls') + KK;
tc_real = 53;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
t_0(1) = 25 + KK;
for i = 2:length(ti)
    t_f = t_e(get_t2(ti(i),v,x));
    t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
end
plot(ti,t_0-KK);
hold on
plot(tim,t);
legend('t_c=48s炉温曲线','t_c=43s炉温曲线','t_c=53s炉温曲线','附件的炉温曲线');
xlabel('时间/s')
ylabel('温度/\circ C')
%% 边界代替t_f的灵敏度分析
t_e = xlsread('z0.xls') + KK;
tc_real = 47.9567;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
t_0(1) = 25 + KK;
for i = 2:length(ti)
    t_f = t_e(get_t2(ti(i),v,x));
    t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
end
plot(ti,t_0-KK);
hold on
t_e = xlsread('z0.xls') + KK;
for i= 1:length(t_e)
    t_e(i) = get_t(x(i)) + KK;
end
tc_real = 43;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
t_0(1) = 25 + KK;
for i = 2:length(ti)
    t_f = t_e(get_t2(ti(i),v,x));
    t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
end
plot(ti,t_0-KK);
hold on
plot(tim,t);
legend('正常拟合的炉温曲线 ','边界温度代替中心温度的炉温曲线','附件的炉温曲线');
xlabel('时间/s')
ylabel('温度/\circ C')

%% 边界温度分布函数
function temp = get_t(x)
to = 25;
t1 = 175;
t2 = 195;
t3 = 235;
t4 = 255;
x = x * 100;
len = 435.5;
if(x<=25)
    temp = to + (t1 - to)/25 * x;
elseif(x<=25 + 30.5*5 + 5*4)
    temp = t1;
elseif(x<=25 + 30.5*5 + 5*5)
    temp = t1 + (t2-t1)/5 * (x - (25 + 30.5*5 + 5*4));
elseif(x<=25 + 30.5*6 + 5*5)
    temp = t2;
elseif(x<=25 + 30.5*6 + 5*6)
    temp = t2 + (t3-t2)/5 * (x - (25 + 30.5*6 + 5*5));
elseif(x<=25 + 30.5*7 + 5*6)
    temp = t3;
elseif(x<=25 + 30.5*7 + 5*7)
    temp = t3 + (t4-t3)/5 * (x - (25 + 30.5*7 + 5*6));
elseif(x<=25 + 30.5*9 + 5*8)
    temp = t4;
else
    temp = (len - x)/(len - (25 + 30.5*9 + 5*8))*(t4 - to) + to;
end
end

%% 时间对应环境温度
function idx = get_t2(tim,v,x)
l = 1;r = length(x);
while(r-l>=5)
    mid = floor((l+r)/2);
    if(tim*v<x(mid))
        r = mid;
    else
        l = mid;
    end
end
for i = l:r
    if(tim*v>=x(i))
        idx = l;
        break;
    end
end
end
