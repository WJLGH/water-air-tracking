function  testfunc()
    % 在工作区声明是需要使用global
    global t 
    % 在函数中无法通过global声明，提升局部变量的作用域
    global q
    t(1);
    q = 23;