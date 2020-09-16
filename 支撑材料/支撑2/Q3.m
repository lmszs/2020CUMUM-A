x0 = [175,195,235,255,80];
lb = [165,185,225,245,65];
ub = [185,205,245,265,100];

options=optimoptions('ga','PopulationSize',100,'MaxGenerations',500);
[x,fval,exitflag,output] = ga(@(x)fun(x),5,[],[],[],[],lb,ub,@(x)nonlcon(x),options)
t1 = x(1);
t2 = x(2);
t3 = x(3);
t4 = x(4);
v = x(5)/100/60;
len = 435.5e-2;
KK = 273.15;
h = 0.5 / 100;
x = 0:h:(len);
% 炉心温度
t_e = zeros(1,length(x));
for i = 1:length(t_e)
    t_e(i) = get_t(t1,t2,t3,t4,x(i)) + KK;
end
% 求解炉温曲线
tc_real = 47.9567;
step = 0.01;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
t_0(1) = 25 + KK;
for i = 2:length(ti)
    t_f = t_e(get_t2(ti(i),v,x));
    t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
end
t_0 = t_0 - KK;
plot(ti,t_0);
check(t_0,ti)


% fun(ub)
%% 
function area = fun(X)
t1 = X(1);
t2 = X(2);
t3 = X(3);
t4 = X(4);
v = X(5)/100/60;
len = 435.5e-2;
KK = 273.15;
h = 0.5 / 100;
x = 0:h:(len);
% 炉心温度
t_e = zeros(1,length(x));
for i = 1:length(t_e)
    t_e(i) = get_t(t1,t2,t3,t4,x(i)) + KK;
end
% 求解炉温曲线
tc_real = 47.9567;
step = 0.01;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
t_0(1) = 25 + KK;
for i = 2:length(ti)
    t_f = t_e(get_t2(ti(i),v,x));
    t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
end
t_0 = t_0 - KK;
[~,idx] = max(t_0);
area = 0;
for i = 1:idx
    if(t_0(i)>217)
        area = area + (t_0(i)-217) * step;
    end
end
end
%% a
function [c,ceq] = nonlcon(X)
t1 = X(1);
t2 = X(2);
t3 = X(3);
t4 = X(4);
v = X(5)/100/60;
len = 435.5e-2;
KK = 273.15;
h = 0.5 / 100;
x = 0:h:(len);
% 炉心温度
t_e = zeros(1,length(x));
for i = 1:length(t_e)
    t_e(i) = get_t(t1,t2,t3,t4,x(i)) + KK;
end
% 求解炉温曲线
tc_real = 47.9567;
step = 0.01;
ti = 0:step:(len/v);     % 时间
t_0 = zeros(1,length(ti));
t_0(1) = 25 + KK;
for i = 2:length(ti)
    t_f = t_e(get_t2(ti(i),v,x));
    t_0(i) = t_0(i-1) + step * (t_0(i-1) - t_f) / (-tc_real);
end
ch = check(t_0-KK,ti);
c = [ch(1) - 3;
     - 3 - ch(2);
     ch(3) - 120;
     60 - ch(3);
     ch(4) - 90;
     40 - ch(4);
     ch(5) - 250;
     240 - ch(5)];
ceq = [];
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
%% 边界温度分布函数
function temp = get_t(t1,t2,t3,t4,x)
to = 25;
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
