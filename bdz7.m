function bdz7(X,Y)
p=0;
syms x y a b c d;
disp('аппроксимация многочленом 1 степени')
for i=1:1:length(X)
p=p+(a*X(i)+b-Y(i))^2;
end
disp('полином:')
p
pa=diff(p,a);disp(['dp/da=']),disp([pa])
pb=diff(p,b);disp(['dp/db=']),disp([pb])
A=[diff(pa,a) diff(pa,b) ;diff(pb,a) diff(pb,b)];
B=-[pa-diff(pa,a)*a-diff(pa,b)*b;pb-diff(pb,a)*a-diff(pb,b)*b;];
XX=A\B;
br=sym2poly(XX(2));
ar=sym2poly(XX(1));
disp(['a = ',num2str(ar), '; b=', num2str(br),';'])
resultl=((subs(subs((a*x+b-y),a,ar),b,br)));
disp('___________________________________________________')
yl=solve(resultl==0,y);
disp('P1(x)='),disp(yl)
disp('___________________________________________________')
q=0:(length(X)-1);
yli=sym2poly(sum(subs(yl,X).*x.^q));
ylsq=Y-fliplr(yli);
sqotln=sqrt((1/length(X))*sum(ylsq.^2));
disp(['среднеквадратическое отклонение = ',num2str(sqotln)])
disp('___________________________________________________')

disp('аппроксимация многочленом 2 степени')
p=0;
for i=1:1:length(X)
p=p+(a*(X(i))^2+b*(X(i))+c-Y(i))^2;
end
disp('полином:')
p
pa=diff(p,a);disp(['dp/da=']),disp([pa])
pb=diff(p,b);disp(['dp/db=']),disp([pb])
pc=diff(p,c);disp(['dp/dc=']),disp([pc])
A=[diff(pa,a) diff(pa,b) diff(pa,c);diff(pb,a) diff(pb,b) diff(pb,c);diff(pc,a) diff(pc,b) diff(pc,c)];
B=-[pa-diff(pa,a)*a-diff(pa,b)*b-diff(pa,c)*c;pb-diff(pb,a)*a-diff(pb,b)*b-diff(pb,c)*c;pc-diff(pc,a)*a-diff(pc,b)*b-diff(pc,c)*c];
XX=A\B;
csq=sym2poly(XX(3));
bsq=sym2poly(XX(2));
asq=sym2poly(XX(1));
disp(['a = ',num2str(asq), '; b=', num2str(bsq),'; c=',num2str(csq),';'])
resultsq=(subs(subs(subs((a*x^2+b*x+c-y),a,asq),b,bsq),c,csq));
disp('___________________________________________________')
ysq=solve(resultsq==0,y);
disp('P2(x)='),disp(ysq)
disp('___________________________________________________')
ysqi=sym2poly(sum(subs(ysq,X).*x.^q));
ysqsq=Y-fliplr(ysqi);
sqotsqn=sqrt((1/length(X))*sum(ysqsq.^2));
disp(['среднеквадратическое отклонение = ',num2str(sqotsqn)])
disp('___________________________________________________')


disp('аппроксимация многочленом 3 степени')
p=0;
for i=1:1:length(X)
p=p+(a*(X(i))^3+b*(X(i))^2+c*(X(i))+d-Y(i))^2;
end
disp('полином:')
p
pa=diff(p,a); disp(['dp/da=']),disp([pa])
pb=diff(p,b);disp(['dp/db=']),disp([pb])
pc=diff(p,c);disp(['dp/dc=']),disp([pc])
pd=diff(p,d);disp(['dp/dd=']),disp([pd])

A=[diff(pa,a) diff(pa,b) diff(pa,c) diff(pa,d);diff(pb,a) diff(pb,b) diff(pb,c) diff(pb,d);diff(pc,a) diff(pc,b) diff(pc,c) diff(pc,d); diff(pd,a) diff(pd,b) diff(pd,c) diff(pd,d)];
B=-[pa-diff(pa,a)*a-diff(pa,b)*b-diff(pa,c)*c-diff(pa,d)*d;pb-diff(pb,a)*a-diff(pb,b)*b-diff(pb,c)*c-diff(pb,d)*d;pc-diff(pc,a)*a-diff(pc,b)*b-diff(pc,c)*c-diff(pc,d)*d;pd-diff(pd,a)*a-diff(pd,b)*b-diff(pd,c)*c-diff(pd,d)*d];
XX=A\B
dcube=sym2poly(XX(4));
ccube=sym2poly(XX(3));
bcube=sym2poly(XX(2));
acube=sym2poly(XX(1));

disp(['a = ',num2str(acube), '; b = ', num2str(bcube),'; c = ',num2str(ccube),'; d = ',num2str(dcube)])
resultsq=subs(subs(subs(subs((a*x^3+b*x^2+c*x+d-y),a,acube),b,bcube),c,ccube),d,dcube);
disp('___________________________________________________')
ycube=solve(resultsq==0,y);
disp('P3(x)='),disp(ycube)
disp('___________________________________________________')
ycubei=sym2poly(sum(subs(ycube,X).*x.^q));
ycubesq=Y-fliplr(ycubei);
sqotcuben=sqrt((1/length(X))*sum(ycubesq.^2));
disp(['среднеквадратическое отклонение = ',num2str(sqotcuben)])
disp('___________________________________________________')

sqots=[sqotln sqotsqn sqotcuben];
[sqotval,I]=min(sqots);
sqotname=["линейный";"квадратный";"кубический"];
disp(['таким образом минимальное с.к.о delta~=',num2str(sqotval)])
disp(['соответствующий ему тип зависимости - ',sqotname(I)])

end

