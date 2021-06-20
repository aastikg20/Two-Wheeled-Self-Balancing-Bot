m_c=1.5;
m_p=0.5;
g=9.82;
L=1;
d1=0.01;
d2=0.01;
A=[0,0,1,0;0,0,0,1;0,g*m_p/m_c,-d1/m_c,-d2/(L*m_c);0,g*(m_c+m_p)/(L*m_c),-d1/(L*m_c),-d2*(m_c+m_p)/((L^2)*m_c*m_p)];
B=[0;0;1/m_c;1/(L*m_c)];
C1=[0,1,0,0];
C2=[1,0,0,0];
D=[0];
sys1=ss(A,B,C1,D);
sys2=ss(A,B,C2,D);
P1=pole(sys1);
P2=pole(sys2);
E1=eig(sys1);
E2=eig(sys2);
[Z1,Gain1]=zero(sys1);
[Z2,Gain2]=zero(sys2);
s1=isstable(sys1);
s2=isstable(sys2);
pzmap(sys1);
pzmap(sys2);
o1=rank(obsv(sys1));
o2=rank(obsv(sys2));
c1=rank(obsv(sys1));
c2=rank(obsv(sys2));
tf(sys1);
tf(sys2);
ss2tf(A,B,C1,D);
ss2tf(A,B,C2,D);
zpk(sys1);
zpk(sys2);

ob1=obsv(sys1);
ob2=obsv(sys2);
co1=ctrb(sys1);
co2=ctrb(sys2);

p=[-1 -0.5 -3 -2];
K= place(A,B,p);

Q=diag([2000 5 4 0.00003]);
R=50000;
N=[0;0;0;0];
K=lqr(sys1,Q,R,N);
A1=A-B*K;
sys_lqr=ss(A1,B,C2,D);
P=pole(sys_lqr);
rlocus(sys_lqr)
S=isstable(sys_lqr)