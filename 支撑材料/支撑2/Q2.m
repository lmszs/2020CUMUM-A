v = (65:100)/60/100;
len = 435.5e-2;
KK = 273.15;

ll = 15e-2;    % 宽度
h = 0.5 / 100;
k = 0.1 / 100;
x = 0:h:(len);
y = 0:k:(ll);
%% 炉心温度
interation = 50000;
n = size(x,2);
m = size(y,2);
z_p = ones(n,m) * 0;
z_c = ones(n,m) * 0;
A = [h^2,h^2,k^2,k^2];
Eps = 1e-3;
w = 1.7;
for int = 1:interation
    for j = 1:m
        z_c(1,j) = 25 + KK;
        z_c(n,j) = 25 + KK;
    end
    for i = 1:n
        z_c(i,m) = get_t(x(i)) + KK;
        z_c(i,1) = get_t(x(i)) + KK;
    end
    for i = 2:n-1
        for j = 2:m-1
            z_c(i,j) = (1-w)*z_p(i,j) + w*(dot(A,[z_c(i,j-1),z_p(i,j+1),z_c(i-1,j),z_p(i+1,j)]))/2/(h^2+k^2);
        end
    end
    err = abs(z_p - z_c);
    if(all(err<Eps))
        break;
    end
    z_p = z_c;
end
z_p = z_p - KK;
tt = ll/2 / k;
xlswrite('z2.xls',z_p(:,tt));
%% 大步长画图
tc_real = 47.9567;
step = 0.01;
t_e = xlsread('z2.xls') + KK;
ch = zeros(5,length(v));
for ttt = 1:length(v)
    ti = 0:step:(len/v(ttt));     % 时间
    t_0 = zeros(1,length(ti));
    t_0(1) = 25 + KK;
    for i = 2:length(ti)
        t_f = t_e(get_t2(ti(i),v(ttt),x));
        t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
    end
    ch(:,ttt) = check(t_0-KK,ti);
end
plot(v*100*60,ch(1,:)*10);
hold on
plot(v*100*60,ch(2,:)*10);
hold on
plot(v*100*60,ch(3,:));
hold on
plot(v*100*60,ch(4,:));
hold on
plot(v*100*60,ch(5,:));
legend('温度上升斜率最高值*10','温度下降斜率最低值*10','上升过程中在150\circC-190\circC的时间','温度大于217?C的时间','峰值温度');
xlabel('速度/(cm/min)')
ylabel('制程界限')
%% 搜索最大可行v
v = (70:0.01:75)/60/100;
tc_real = 47.9567;
step = 0.01;
t_e = xlsread('z2.xls') + KK;
ch = zeros(5,length(v));
for ttt = 1:length(v)
    ti = 0:step:(len/v(ttt));     % 时间
    t_0 = zeros(1,length(ti));
    t_0(1) = 25 + KK;
    for i = 2:length(ti)
        t_f = t_e(get_t2(ti(i),v(ttt),x));
        t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
    end
    ch(:,ttt) = check(t_0-KK,ti);
end
res = 0;
for i = length(v):-1:1
    if(ch(1,i)<=3 && ch(2,i)>=-3 && ch(3,i)>=60 && ch(3,i)<=120 &&...
        ch(4,i)>=40 && ch(4,i)<=90 && ch(5,i)>=240 && ch(5,i)<=250)
        res = v(i);
        break
    end
end
res*60*100
%% 最优v=72.25画图
v = 72.25/60/100;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
t_0(1) = 25 + KK;
for i = 2:length(ti)
    t_f = t_e(get_t2(ti(i),v,x));
    t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
end
plot(ti,t_0-KK);
xlabel('时间/s')
ylabel('温度/\circ C')
ch = check(t_0-KK,ti)


%% 炉温曲线检验
function res = check(t, ti)
delta = ti(2) - ti(1);
res = [0,1e9,0,0,0];
for i = 2:length(t)
    if(i>1)
        res(1) = max(res(1),(t(i)-t(i-1))/delta);
        res(2) = min(res(2),(t(i)-t(i-1))/delta);
        if(t(i)-t(i-1)>=0 && t(i)<=190 && t(i)>=150)
            res(3) = res(3) + delta;
        end
    end
    if(t(i)>217)
        res(4) = res(4) + delta;
    end
end
res(5) = max(t);
end
%% 边界温度分布函数
function temp = get_t(x)
to = 25;
t1 = 182;
t2 = 203;
t3 = 237;
t4 = 254;
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

%% 时间对应温度
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
