%% TASK 1

m_c=1.5;
m_p=0.5;
g=9.82;
L=1;
d1=0.01;
d2=0.01;
A=[0,0,1,0;0,0,0,1;0,g*m_p/m_c,-d1/m_c,-d2/(L*m_c);0,g*(m_c+m_p)/(L*m_c),-d1/(L*m_c),-d2*(m_c+m_p)/((L^2)*m_c*m_p)];
B=[0;0;1/m_c;1/(L*m_c)];
c=[0,1,0,0];
D=[0];
sys=ss(A,B,c,D);
P=pole(sys);
E=eig(sys);
[Z,G]=zero(sys);
pzmap(sys);
tf(sys);
ss2tf(A,B,c,D);
zeros=[2.933381967792379e-17;0];
poles=[-3.6327,3.6043,-0.0050,0];
gain=0.6667;
sys=zpk(zeros,poles,gain);  

%% TASK 2

[A,B,c,D]=ssdata(sys);
rlocus(sys);
p=[-3.632666683606961;-3.604333347091377;-0.004999996817749;0];
K=place(A,B,p);
sys_c1=ss(A-B*K,B,c,D);
sisotool(sys_c1);
