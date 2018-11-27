clc;
clear;
close all;

x = rand(100,1);
y = rand(100,1);
data0 = [x y];%% data
[a b]=size(data0);
d=zeros(a,1);
f1=1;
f2=1;

for i = 1:a
    if y(i) < 0.2 * sin(10 * x(i)) + 0.3 || (y(i) - 0.8)^2 + (x(i) - 0.5)^2 < 0.15^2
        d(i,1) = 1;
        data1(f1,:)=[x(i) y(i)];%% data 1
        f1=f1+1;
    else
        d(i,1) = -1;
        data2(f2,:)=[x(i) y(i)];%% data 2
        f2=f2+1;
    end
end
for i = 1:a
    if d(i,1) == 1
        scatter(data0(i,1),data0(i,2),'g*');
        axis([0,1,0,1]);
        hold on;
    else
      scatter(data0(i,1),data0(i,2),'r+');
      axis([0,1,0,1]);
      hold on;
    end
end

%%boundary
n = 0:0.001:1;
y1 = 0.2 * sin(10 * n) + 0.3;
plot(n,y1,'k');
axis([0,1,0,1]);
hold on;
rectangle('Position',[0.35,0.65,0.3,0.3],'Curvature',[1,1]);
hold on

%%% k_means
%k_means(data1,10) ;  
   [u1 re]=KMeans(data1,2); 
   [u2 re]=KMeans(data2,2); 
 
% % ??????????

hold on;

 
scatter(u1(:,1),u1(:,2),'ro'); 
hold on;
scatter(u2(:,1),u2(:,2),'bo');   
hold on;

p=0;q=0;
for i=1:2
    if isnan(u1(i,1))~=1
       p=p+1;
       c1(p,:)=u1(i,:);
    end
    if isnan(u2(i,1))~=1
       q=q+1;
       c2(q,:)=u2(i,:);
    end
 
end


%% g(x)
c=[c1;c2];
w=rand(p+q,1)-0.5;
theta=rand(100,1)-0.5;


eta=0.05;
gx=zeros(100,1);
epoch=1;
errors=zeros(epoch,1);
lamda=0.000001;
while epoch<200
    errors(epoch)=0;
for ii=1:100
    for i=1:p+q
    gx(ii)=w(i)*exp(-(norm(data0(ii,:)-c(i,:)))^2*lamda)+gx(ii);
    end
    gx(ii)=gx(ii)+theta(ii);
    if gx(ii)>=0
        gx(ii)=1;
    else
        gx(ii)=-1;
    end
   
    if d(ii)~=gx(ii)
        errors(epoch)=errors(epoch)+1;
    end
end

for ii=1:100
  w=w+eta*exp(-(norm(data0(ii,:)-c(i,:)))^2*lamda)*(d(ii)-gx(ii));
  theta(ii)=theta(ii)+(d(ii)-gx(ii))*eta;
end


epoch=1+epoch;
end

%%PTA
figure;
xx=1:epoch-1;
plot(xx,errors);

%%contour
figure
scatter(c1(:,1),c1(:,2),'ro'); 
hold on;
scatter(c2(:,1),c2(:,2),'bo'); 
hold on;


%%%%%%%%%%%%
%contour(c1(:,1),c1(:,2),c1);
gridSize = 8;
u = linspace(0, 1, gridSize);
v = linspace(0, 1, gridSize);
scores1 = zeros(length(u), length(v));
scores2 = zeros(length(u), length(v));
pp = zeros(length(u), length(v));

 for iii=1:100
            if d(iii)==1
                a=ceil(data0(iii,1)*gridSize);
                b=ceil(data0(iii,2)*gridSize);
                scores1(a,b)=scores1(a,b)+1;
            elseif d(iii)==-1
                a=ceil(data0(iii,1)*gridSize);
                b=ceil(data0(iii,2)*gridSize);
                scores2(a,b)=scores2(a,b)+1;
            end
end
for (i = 1 : length(u))
    for j = 1 : length(v)
       
        if (scores1(i, j) == scores2(i, j))
            p(i,j) = 1.5;
        elseif (scores1(i, j) > scores2(i, j))
            p(i, j) = 1;
        else 
            p(i, j) = 2;
        end
    end
        
end
        
    
hold on;

[c, h] = contour(u, v, p',[1.5 1.5]);

hold on;
for i = 1:100
    if d(i,1) == 1
        scatter(data0(i,1),data0(i,2),'g*');
        axis([0,1,0,1]);
        hold on;
    else
      scatter(data0(i,1),data0(i,2),'r+');
      axis([0,1,0,1]);
      hold on;
    end
end


