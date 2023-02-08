pos = [1,0.5];
t = 0:0.1:10;
ints = zeros(length(t),2);
for i = 1:length(t)
    ints(i,:) = get_intensity(pos,t(i));
end

mean(ints)
std(ints)
plot(ints(:,1));
hold on
plot(ints(:,2));
hold off
legend('p1','p2')