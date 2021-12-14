ratio1 = zeros(Nc,1);
ratio = zeros(Nc,1);
for i = 1:Nc
    ratio1(i) = sum(sum(abs(H(:,i)*u(i,:))./abs(y)))/(M*L);
    ratio(i) = sum(sum(abs(G(:,i)*message(i,:))./abs(y)))/(M*L);
end
ratio1 = sort(ratio1);
ratio = sort(ratio);
x=(1:512);
figure;

plot(x,ratio,x,ratio1);
xlabel('用户编号');
ylabel('到达信号占接收信号的比值');
legend('真实比值','仿真比值')