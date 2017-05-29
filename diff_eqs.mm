Func Matrix diff_eqs_nonliner(t,xx,u) // 非線形モデル
Real t;
Matrix xx,u;
{
  Matrix K; // 行列
  Matrix xxdot,z;
  Matrix tmp,tmp2,tmp3;
  Real a,m,M,g,l,f,c,J;

  a=0.73;   // [N/V] uから台車に働く力までのゲイン
  m=0.038;  // [kg] 振子の質量
  M=1.49;   // [kg] 台車・プーリ・モータ系の等価質量
  g=9.8;    // [m/s^2] 重力加速度
  l=0.13;   // [m] 軸から振子の重心までの距離
  f=15.10;  // [kg/s]
  c=2.1e-4; // [kgm^2/s] 軸の摩擦係数
  J=4.5e-4; // [kgm^2] 振子の重心周りの慣性モーメント

  u=[0];

  K=[[ M+m,           m*l*cos(xx(2)) ]
     [ m*l*cos(xx(2)),    J+m*l*l    ]];

  tmp=[[ a*u ]
       [ 0.0 ]];

  tmp2=[[ -m*l*sin(xx(2))*(xx(4)^2) ]
        [      m*g*l*sin(xx(2))    ]];

  tmp3=[[ f*xx(3) ]
        [ c*xx(4) ]];

  z=(K~)*(tmp+tmp2-tmp3);


  xxdot=trans([xx(3) xx(4) z(1) z(2)]); //転置行列

  return xxdot;
}

Func Matrix diff_eqs_liner(t,xx,u) // 線形モデル
Real t;
Matrix xx,u;
{
  Matrix K,_K,A,B;
  Matrix tmp,tmp2,tmp3;
  Matrix xxdot;
  Real a,m,M,g,l,f,c,J;

  a=0.73;   // [N/V] uから台車に働く力までのゲイン
  m=0.038;  // [kg] 振子の質量
  M=1.49;   // [kg] 台車・プーリ・モータ系の等価質量
  g=9.8;    // [m/s^2] 重力加速度
  l=0.13;   // [m] 軸から振子の重心までの距離
  f=15.10;  // [kg/s]
  c=2.1e-4; // [kgm^2/s] 軸の摩擦係数
  J=4.5e-4; // [kgm^2] 振子の重心周りの慣性モーメント

  u=[0];

  K=[[ M+m,      -m*l   ]
     [ -m*l,    J+m*l*l ]];

  /** inv(K) */
  _K=K~;

  tmp=_K*[[ 0,    0    ]
          [ 0,  -m*g*l ]];

  tmp2=_K*[[ -f,   0  ]
           [  0,   -c ]];

  A=[[ Z(2,2) I(2,2) ]
     [ tmp    tmp2   ]];
  // print A;

  tmp3=_K*[[a] [0]];
  B=trans([0 0 tmp3(1) tmp3(2)]);
  // print B;

  xxdot=A*xx+B*u;

  return xxdot;
}

Func void main()
{
  Real t0,t1,h;
  Matrix x0,T,X;

  x0=trans([0 10.0*PI/(180) 0 0]);
  t0=0.0;
  t1=10.0;
  h=0.1;

  // print diff_eqs_nonliner(1,x0,[0]);

  {T,X}=Ode(t0,t1,x0,diff_eqs_liner,"",h);
  mgplot(1,X);
  mgplot_eps(1,"liner_under.eps");
}
