

f = @(theta_1, theta_3) [30*cos(theta_1 + theta_3 + pi/12) + 250*cos(theta_1) + 100*cos(theta_1 + pi/12);
                         30*sin(theta_1 + theta_3 + pi/12) + 250*sin(theta_1) + 100*sin(theta_1 + pi/12)];
                    
x = [350; 30];% 所需的 xy 位置

eqs = @(theta) f(theta(1), theta(2)) - x;
theta0 = [0; 0];
theta_sol = fsolve(eqs, theta0);

theta1_sol = theta_sol(1);
theta3_sol = theta_sol(2);