function [M,fLagr,difmax,omega,rpoly,rmax] = bdz5( X, Y,forig) %% X - значения аргумента, Y - значениня функции в точках X, 
% forig - функция подвергаемая интерполяции (нужна для вычисления погрешностие
disp('матрица Вардемонда');
A=repmat(X.', 1, length(X)).^repmat([0:length(X)-1], length(X), 1) %объяснение есть в методичке
Y=Y';
disp('строим матрицу А^(-1)')
A=inv(A)
disp('вектор коэфициентов для полинома')
M=(A*Y) %решаем уравнение относительно коэфффффициентов

syms x;
for i=1:length(M)
    pwr(i)=(length(M)-i);%задаю степени полинома
end
syms x;
M=flipud(M)';
disp('коэфициента полинома для вывода по степеням. Например: M(1)*x^3')
disp(M)
disp('полином лагранжа ');
fLagr=sum(M.*(x.^pwr),2) %строю полином
xx=X(1):0.1:X(length(X)); %сетка для вычислений
disp('(n+1) производная интерполируемой функции (нужна для оценки погрешности)');
difforig=diff(forig,length(X))
[difmax,I]=max(subs(diff(forig,length(X)),xx));
disp('максимум (n+1)й производной в числителе равен ')
difmax=sym2poly(difmax)
disp('и достигается при x равном/приблизительно равном')
disp(xx(I))
omega=1;
for i=1:1:length(X)
    omega=omega*(x-X(i));
end
disp('полином омега(n+1)')
omega
rpoly=(difmax/factorial(length(X)))*omega;
[rmax,I]=max(subs(rpoly,xx));
disp('максимум погрешности');
rmax=sym2poly(rmax)
disp('приблизительное значение х в котором погрешность максимальна, указывать не нужно')
disp(xx(I))
end 
