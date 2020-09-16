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
%% 炉心温度
temp = zeros(1,length(x));
for i = 1:length(x)
    temp(i) = get_t(x(i));
end
plot(x*100,temp)
hold on
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
    err = sqrt(sum((z_p - z_c).^2)/length(z_c));
    if(err<Eps)
        break;
    end
    z_p = z_c;
end
z_p = z_p - KK;
% [x,y]=meshgrid(x,y);
% mesh(x',y',z_p);
tt = ll/2 / k;
xlswrite('z0.xls',z_p(:,tt));
plot(x*100,z_p(:,tt))
legend('边界线温度','中心线温度')
xlabel('位置/cm')
ylabel('温度/\circC')
%% 大范围搜索t_c
step = 0.01;
t_e = xlsread('z0.xls') + KK;
hold on
ti = 0:step:(len/v);     % 时间
t_c = 1:1:100;
t_0 = zeros(1,length(ti));
err = zeros(1,length(t_c));
for kkkk = 1:length(t_c)
    t_0(1) = 25 + KK;
    for i = 2:length(ti)
        t_f = t_e(get_t2(ti(i),v,x));
        t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-t_c(kkkk));
    end
    ccc = 1;
    err(kkkk) = 0;
    for i = 1:length(ti)
        if(ccc>length(t))
            break;
        end
        if(abs(ti(i) - tim(ccc))<step/2)
            err(kkkk) = err(kkkk) + (t_0(i)-KK - t(ccc))^2;
            ccc = ccc + 1;
        end
    end
    err(kkkk) = err(kkkk) / length(t);
    err(kkkk) = sqrt(err(kkkk));
end
plot(t_c,err)
xlabel('时间常数t_c/s');
ylabel('平均误差/\circC');
%% 三分法
l = 0.01;r = 100;
step = 0.01;
t_e = xlsread('z0.xls') + KK;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
while(r-l>=1e-6)
    m1 = (2*l+r)/3;
    m2 = (l+2*r)/3;
    % e1
    t_0(1) = 25 + KK;
    for i = 2:length(ti)
        t_f = t_e(get_t2(ti(i),v,x));
        t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-m1);
    end
    ccc = 1;
    e1 = 0;
    for i = 1:length(ti)
        if(ccc>length(t))
            break;
        end
        if(abs(ti(i) - tim(ccc))<step/2)
            e1 = e1 + (t_0(i)-KK - t(ccc))^2;
            ccc = ccc + 1;
        end
    end
    e1 = e1 / length(t);
    e1 = sqrt(e1);
    re = check(t_0,ti);
    if(re(5)-KK<240) 
        e1 = 1e7;
    end
    % e2
    t_0(1) = 25 + KK;
    for i = 2:length(ti)
        t_f = t_e(get_t2(ti(i),v,x));
        t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-m2);
    end
    ccc = 1;
    e2 = 0;
    for i = 1:length(ti)
        if(ccc>length(t))
            break;
        end
        if(abs(ti(i) - tim(ccc))<step/2)
            e2 = e2 + (t_0(i)-KK - t(ccc))^2;
            ccc = ccc + 1;
        end
    end
    e2 = e2 / length(t);
    e2 = sqrt(e2);
    re = check(t_0,ti);
    if(re(5)-KK<240) 
        e2 = 1e7;
    end
    if(e1>e2)
        l = m1;
    else
        r = m2;
    end
end
l
r
e1
e2
%% 最优t_c = 47.9567画图(err = 9.55)
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
plot(tim,t);
legend('t_c=47.9567s炉温曲线','附件的炉温曲线');
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