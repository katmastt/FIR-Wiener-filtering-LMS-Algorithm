load('1.mat')

p1=10;
autr_x_all=xcorr(x,p1-1);
autr_x=autr_x_all(p1:2*p1-1);
R_x_inv = inv(toeplitz(autr_x));
autr_d_all=xcorr(d,p1-1);
autr_d=transpose(autr_d_all(p1:2*p1-1));
w=transpose(R_x_inv*autr_d);
d_exp_m1=zeros(1,300);
for i=2:301
    for c=1:p1
       if i-c>=2
            d_exp_m1(i-1)=d_exp_m1(i-1)+w(c)*x(i-1-c);
       elseif i-c==1
            d_exp_m1(i-1)=d_exp_m1(i-1)+w(c)*x(1);
       else
            d_exp_m1(i-1)=d_exp_m1(i-1)+w(c)*x(-i+c+1);
       end
    end
end

plot(d,'k-','LineWidth',2); grid on; hold on;
plot(d_exp_m1,'b-','LineWidth',2);
xlabel('Sample index');
ylabel('Signal Values'); 
axis([0 300 -2 8]);
legend('Desired signal d[n]','Filtered signal with p=10', 'Location','NorthEast')
title('Filtering via Wiener-Hopf method');
figure();

p2=40;
autr_x_all=xcorr(x,p2-1);
autr_x=autr_x_all(p2:2*p2-1);
R_x_inv = inv(toeplitz(autr_x));
autr_d_all=xcorr(d,p2-1);
autr_d=transpose(autr_d_all(p2:2*p2-1));
w=transpose(R_x_inv*autr_d);
d_exp_m1=zeros(1,300);
for i=2:301
    for c=1:p2
       if i-c>=2
            d_exp_m1(i-1)=d_exp_m1(i-1)+w(c)*x(i-1-c);
       elseif i-c==1
            d_exp_m1(i-1)=d_exp_m1(i-1)+w(c)*x(1);
       else
            d_exp_m1(i-1)=d_exp_m1(i-1)+w(c)*x(-i+c+1);
       end
    end
end
plot(d,'k-','LineWidth',2); grid on; hold on;
plot(d_exp_m1,'b-','LineWidth',2);
xlabel('Sample index');
ylabel('Signal Values'); 
axis([0 300 -2 8]);
legend('Desired signal d[n]','Filtered signal with p=40', 'Location','NorthEast')
title('Filtering via Wiener-Hopf method');
figure();

p3=80;
autr_x_all=xcorr(x,p3-1);
autr_x=autr_x_all(p3:2*p3-1);
R_x_inv = inv(toeplitz(autr_x));
autr_d_all=xcorr(d,p3-1);
autr_d=transpose(autr_d_all(p3:2*p3-1));
w=transpose(R_x_inv*autr_d);
d_exp_m1=zeros(1,300);
for i=2:301
    for c=1:p3
       if i-c>=2
            d_exp_m1(i-1)=d_exp_m1(i-1)+w(c)*x(i-1-c);
       elseif i-c==1
            d_exp_m1(i-1)=d_exp_m1(i-1)+w(c)*x(1);
       else
            d_exp_m1(i-1)=d_exp_m1(i-1)+w(c)*x(-i+c+1);
       end
    end
end
plot(d,'k-','LineWidth',2); grid on; hold on;
plot(d_exp_m1,'b-','LineWidth',2);
xlabel('Sample index');
ylabel('Signal Values'); 
axis([0 300 -2 8]);
legend('Desired signal d[n]','Filtered signal with p=80', 'Location','NorthEast')
title('Filtering via Wiener-Hopf method');


clear
load('2.mat')

m1=0.012;
m2=0.05;
p3=2;
w_n_m1=zeros(300,p3);
w_n_m2=zeros(300,p3);
d_exp_m1=zeros(1,300);
e_m1=zeros(1,300);
d_exp_m2=zeros(1,300);
e_m2=zeros(1,300);
x=zeros(1,p3);
for i=1:300
        p=1;
        index=1;
        while i-p>=0 && index<=p3
            x(1,index)=u(i-p+1);
            index=index+1;
            p=p+1;
        end
            d_exp_m1(i)=w_n_m1(i,:)*transpose(x(1,:));
            e_m1(i)=y(i)-d_exp_m1(i);
            d_exp_m2(i)=w_n_m2(i,:)*transpose(x(1,:));
            e_m2(i)=y(i)-d_exp_m2(i);
            w_n_m1(i+1,:)=w_n_m1(i,:)+m1*e_m1(i)*x(1,:); 
            w_n_m2(i+1,:)=w_n_m2(i,:)+m2*e_m2(i)*x(1,:); 
   
end
figure();
plot(w_n_m1(:,1),'k-','LineWidth',2); grid on; hold on;
plot(w_n_m2(:,1),'b-','LineWidth',2); grid on; hold on;
plot(w_n_m1(:,2),'k-','LineWidth',2,'LineStyle','--'); grid on; hold on;
plot(w_n_m2(:,2),'b-','LineWidth',2,'LineStyle','--'); grid on; hold on;

title('Adaptive Filtering via LMS with p=2');
xlabel('Iteration n');
ylabel('LMS-Based FIR Filter Coefficients');
axis([0 300 -4 4]);
legend('W_n[0] for m=0.012','W_n[0] for m=0.05','W_n[1] for m=0.012','W_n[1] for m=0.05', 'Location','NorthEast')


p2=3;
w_n_m1=zeros(300,p2);
w_n_m2=zeros(300,p2);
d_exp_m1=zeros(1,300);
e_m1=zeros(1,300);
d_exp_m2=zeros(1,300);
e_m2=zeros(1,300);
x=zeros(1,p2);
for i=1:300
        p=1;
        index=1;
        while i-p>=0 && index<=p2
            x(1,index)=u(i-p+1);
            index=index+1;
            p=p+1;
        end
            d_exp_m1(i)=w_n_m1(i,:)*transpose(x(1,:));
            e_m1(i)=y(i)-d_exp_m1(i);
            d_exp_m2(i)=w_n_m2(i,:)*transpose(x(1,:));
            e_m2(i)=y(i)-d_exp_m2(i);
            w_n_m1(i+1,:)=w_n_m1(i,:)+m1*e_m1(i)*x(1,:); 
            w_n_m2(i+1,:)=w_n_m2(i,:)+m2*e_m2(i)*x(1,:); 
   
end
figure();
plot(w_n_m1(:,1),'k-','LineWidth',2); grid on; hold on;
plot(w_n_m2(:,1),'b-','LineWidth',2); grid on; hold on;
plot(w_n_m1(:,2),'k-','LineWidth',2,'LineStyle','--'); grid on; hold on;
plot(w_n_m2(:,2),'b-','LineWidth',2,'LineStyle','--'); grid on; hold on;
plot(w_n_m1(:,3),'k-','LineWidth',2,'LineStyle','-.'); grid on; hold on;
plot(w_n_m2(:,3),'b-','LineWidth',2,'LineStyle','-.'); grid on; hold on;
title('Adaptive Filtering via LMS with p=3');
xlabel('Iteration n');
ylabel('LMS-Based FIR Filter Coefficients');
axis([0 300 -3 4]);
legend('W_n[0] for m=0.012','W_n[0] for m=0.05','W_n[1] for m=0.012','W_n[1] for m=0.05','W_n[2] for m=0.012','W_n[2] for m=0.05', 'Location','NorthEast')



p1=4;
w_n_m1=zeros(300,p1);
w_n_m2=zeros(300,p1);
d_exp_m1=zeros(1,300);
e_m1=zeros(1,300);
d_exp_m2=zeros(1,300);
e_m2=zeros(1,300);
x=zeros(1,p1);
for i=1:300
        p=1;
        index=1;
        while i-p>=0 && index<=p1
            x(1,index)=u(i-p+1);
            index=index+1;
            p=p+1;
        end
            d_exp_m1(i)=w_n_m1(i,:)*transpose(x(1,:));
            e_m1(i)=y(i)-d_exp_m1(i);
            d_exp_m2(i)=w_n_m2(i,:)*transpose(x(1,:));
            e_m2(i)=y(i)-d_exp_m2(i);
            w_n_m1(i+1,:)=w_n_m1(i,:)+m1*e_m1(i)*x(1,:); 
            w_n_m2(i+1,:)=w_n_m2(i,:)+m2*e_m2(i)*x(1,:); 
   
end
figure();
plot(w_n_m1(:,1),'k-','LineWidth',2); grid on; hold on;
plot(w_n_m2(:,1),'b-','LineWidth',2); grid on; hold on;
plot(w_n_m1(:,2),'k-','LineWidth',2,'LineStyle','--'); grid on; hold on;
plot(w_n_m2(:,2),'b-','LineWidth',2,'LineStyle','--'); grid on; hold on;
plot(w_n_m1(:,3),'k-','LineWidth',2,'LineStyle','-.'); grid on; hold on;
plot(w_n_m2(:,3),'b-','LineWidth',2,'LineStyle','-.'); grid on; hold on;
plot(w_n_m1(:,4),'k-','LineWidth',2,'LineStyle',':'); grid on; hold on;
plot(w_n_m2(:,4),'b-','LineWidth',2,'LineStyle',':'); grid on; hold on;
title('Adaptive Filtering via LMS with p=4');
xlabel('Iteration n');
ylabel('LMS-Based FIR Filter Coefficients');
axis([0 300 -3 4]);
legend('W_n[0] for m=0.012','W_n[0] for m=0.05','W_n[1] for m=0.012','W_n[1] for m=0.05','W_n[2] for m=0.012','W_n[2] for m=0.05','W_n[3] for m=0.012','W_n[3] for m=0.05', 'Location','NorthEast')

