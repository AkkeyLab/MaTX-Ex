Func Matrix diff_eqs_nonliner(t,xx,u) // 非線形モデル
Real t;
Matrix xx,u;
{
  Matrix K; // 行列
  Matrix xxdot, z; // 計算結果格納と途中結果格納用
  Matrix _t1, _t2, _t3; // 途中計算用
  Real a, m, M, g, l, f, c, J;

  a = 0.73;   // [N/V] uから台車に働く力までのゲイン
  m = 0.038;  // [kg] 振子の質量
  M = 1.49;   // [kg] 台車・プーリ・モータ系の等価質量
  g = 9.8;    // [m/s^2] 重力加速度
  l = 0.13;   // [m] 軸から振子の重心までの距離
  f = 15.10;  // [kg/s]
  c = 2.1e-4; // [kgm^2/s] 軸の摩擦係数
  J = 4.5e-4; // [kgm^2] 振子の重心周りの慣性モーメント

  u = [0];    // パワーアンプへの入力電圧（バグ対策）

/**
  x = [r θ r' θ']^t  ->  x' = [r' θ' r'' θ'']^t

  # Process
  r''とθ''を求める必要がある（それ以外はxx(i)として参照して利用）
  -> K^-1[-fr'+mlsinθ*θ'^2+au mglsinθ-cθ']^t
  -> (K = [M+m mlcosθ],[mlcosθ J+ml^2])
  [-fr'+mlsinθ*θ'^2+au mglsinθ-cθ']^t
  = [mlsinθ*θ'^2 mlsinθ*g] - [fr' cθ'] + [au 0]
  = [mlsin*xx(2)*xx(4)^2 mlsin*xx(2)*g] - [f*xx(3) c*xx(4)] + [au 0]
*/

  K = [[ M + m,              m * l * cos(xx(2)) ]
      [  m * l * cos(xx(2)), J + m * l * l      ]];

  _t1 = [[ a * u ]
         [ 0.0   ]];

  _t2 = [[ -m * l * sin(xx(2)) * (xx(4)^2) ] // 資料とは制御が多少異なるため調整
         [  m * l * sin(xx(2)) * g         ]];

  _t3 = [[ f * xx(3) ]
         [ c * xx(4) ]];

  z = (K~) * (_t1 + _t2 - _t3);


  xxdot = trans([xx(3) xx(4) z(1) z(2)]); //転置行列

  return xxdot;
}

Func Matrix diff_eqs_liner(t,xx,u) // 線形モデル
Real t;
Matrix xx, u;
{
  Matrix K, _K, A, B;
  Matrix _t1, _t2, _t3;
  Matrix xxdot;
  Real a, m, M, g, l, f, c, J;

  a = 0.73;   // [N/V] uから台車に働く力までのゲイン
  m = 0.038;  // [kg] 振子の質量
  M = 1.49;   // [kg] 台車・プーリ・モータ系の等価質量
  g = 9.8;    // [m/s^2] 重力加速度
  l = 0.13;   // [m] 軸から振子の重心までの距離
  f = 15.10;  // [kg/s]
  c = 2.1e-4; // [kgm^2/s] 軸の摩擦係数
  J = 4.5e-4; // [kgm^2] 振子の重心周りの慣性モーメント

  u = [0];    // パワーアンプへの入力電圧（バグ対策）

  K = [[  M + m,    -m * l     ]
       [ -m * l, J + m * l * l ]];

  _K = K~; // 別表記：inv(K)

  _t1 = _K * [[ 0, 0          ]
              [ 0, -m * g * l ]];

  _t2 = _K * [[ -f, 0  ]
              [  0, -c ]];

  A=[[ Z(2,2) I(2,2) ]
     [ _t1    _t2    ]];

  _t3 = _K*[[a] [0]];

  B = trans([0 0 _t3(1) _t3(2)]);

  xxdot = A * xx + B * u;

  return xxdot;
}

Func void main()
{
  Real t0, t1, h;
  Matrix x0, T, X, Y;

  x0 = trans([0 10.0 * PI / (180) 0 0]);
  t0 = 0.0;
  t1 = 10.0;
  h  = 0.1;

  {T, Y} = Ode(t0, t1, x0, diff_eqs_nonliner, "", h);
  mgplot(1, Y);
  mgplot_eps(1, "diff_eqs_nonliner.eps");

  {T, X} = Ode(t0, t1, x0, diff_eqs_liner, "", h);
  mgplot(1, X);
  mgplot_eps(1, "diff_eqs_liner.eps");
}
