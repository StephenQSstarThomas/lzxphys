# 群元素 Li₂ 主公式与 E 型精确几何（2026-09-23）

本稿执行用户重新指定的路线：**Li₂ 解析式与纯符号化简优先，直接几何为后续方法**。
不使用陪集转移，不进行数值拟合或数值归约，也不把另一个同类代表冒充原几何相位。
09-22 的 transfer 设计已明确标记为被否决；旧计数代表仅保留作历史材料，不参与本稿计算。

范围为 A=循环群（沿用）、D=双重二面体群（沿用）、E=2T/2O/2I。
本稿新增一般位置的纯符号主程序、三个 E 群的完整非退化几何分类、
2T 的十二类实际有理体积，以及覆盖退化 E 输入的原 E/F/T 精确解析接口。
必须区分“完整解析表达”与“全部特殊函数已被消去”：后者对 2O/2I 并未成立，
而“所有体积都是 pi² 的有理倍数”已被本稿的严格反例排除。

## 1. 六个迹、取向和适用范围

采用 g=x0 I+i(x1 sigma1+x2 sigma2+x3 sigma3)，四元数乘法为负叉乘约定。
从三个输入构造

\[
(p_0,p_1,p_2,p_3)=(I,g_1,g_1g_2,g_1g_2g_3),\qquad V=(q(p_0),q(p_1),q(p_2),q(p_3)).
\]

本稿严格使用用户给出的边序，不改动旧数值 API 的历史边序：

\[
\begin{aligned}
\ell_1&=\arccos\tfrac12\operatorname{Tr}g_1,&
\ell_2&=\arccos\tfrac12\operatorname{Tr}(g_1g_2),&
\ell_3&=\arccos\tfrac12\operatorname{Tr}(g_1g_2g_3),\\
\ell_4&=\arccos\tfrac12\operatorname{Tr}g_3,&
\ell_5&=\arccos\tfrac12\operatorname{Tr}(g_2g_3),&
\ell_6&=\arccos\tfrac12\operatorname{Tr}g_2.
\end{aligned}
\]

即 (01,02,03,23,13,12)，对边为 j 与 j+3。Gram 矩阵按这些顶点直接取内积。
令 A=g1、B=g1g2、C=g1g2g3，则

\[
D=\det V=\tfrac14\operatorname{Tr}[A(BC-CB)].
\]

这是累计顶点的交换子，**不是**三个增量 g1/g2/g3 的直接交换子。
由泡利乘法 BC-CB=-2i(b×c)·sigma 得到右边=a·(b×c)，即所需行列式。
例如输入 e1、e2、(0,3/5,0,4/5) 时，正确 D=3/5；误用三个增量会得到 4/5。

D≠0 保证四向量独立，因而其正锥不含零；存在共同开半球，唯一凸短测地线四面体成立。
单四面体入口拒绝 D=0，不能用 sgn(0)=0 将全局反常错误置为 1。
三个群元同属 U(1) 时恰好都处在大圆内，所以 A 型限制必须保留其原有退化链/进位构造。

## 2. 完整主支解析式

以下是 [Murakami 定理 1.2](https://arxiv.org/pdf/1011.2584v4) 的群元代入及固定 z 偏导展开。
令 c_j=cos(ell_j)、b_j=c_j+i sqrt(1-c_j²)，
a_j=-b_(j+3)^(-1)，所有角在 (0,pi)，平方根为正根。定义

\[
\begin{aligned}
q_0={}&a_1a_4+a_2a_5+a_3a_6+a_1a_2a_6+a_1a_3a_5+a_2a_3a_4+a_4a_5a_6+\prod_{j=1}^6a_j,\\
q_1={}&-\sum_{j=1}^3(a_j-a_j^{-1})(a_{j+3}-a_{j+3}^{-1}),\qquad q_2=\overline{q_0},\\
\Delta={}&q_1^2-4q_0q_2=16\det G,\qquad
z=\frac{-q_1+\sqrt\Delta}{2q_2}=\frac{-2q_0}{q_1+4\sqrt{\det G}}.
\end{aligned}
\]

令 P=(1245,1346,2356)、N=(456,234,135,126)。串号表示边指标集合。

\[
Z_0=z,\quad (Z_1,Z_2,Z_3)=\left(z\prod_{j\in S}b_j\right)_{S\in P},\quad
(Z_4,Z_5,Z_6,Z_7)=\left(z\prod_{j\in S}b_j\right)_{S\in N}.
\]

注意 N 是**代入边长之后**的支持；不能与原 a 变量中的 123/156/246/345 混淆。
对 j=1,...,6，全部偏导通过下式给出，不留下隐式求导：

\[
\Sigma_j=\sum_{r:j\in P_r}\log(1-Z_r)
-\sum_{r:j\in N_r}\log(1-Z_{r+3}),\qquad
d_j=\frac{\pi-\ell_{j_{\rm opp}}+\operatorname{Im}\Sigma_j}{2},
\]

其中 P_r 的 r=1,2,3，N_r 的 r=1,2,3,4。六行正负支持依次为
(12;67)、(13;57)、(23;56)、(12;45)、(13;46)、(23;47)。

\[
\begin{aligned}
W={}&\frac12\operatorname{Re}\left[\sum_{r=0}^3\operatorname{Li}_2(Z_r)
-\sum_{r=4}^7\operatorname{Li}_2(Z_r)\right]
-\frac12\sum_{j=1}^3(\pi-\ell_j)(\pi-\ell_{j+3})\\
&-\pi\operatorname{Arg}(-q_2)-\sum_{j=1}^6\ell_jd_j-\frac{\pi^2}{2},\\
\omega_k(g_1,g_2,g_3)={}&\exp\left[-\frac{ik}{\pi}\operatorname{sgn}(D)W\right],\qquad k\in\mathbb Z.
\end{aligned}
\]

所有 Li₂、log、Arg 使用主支；W 是体积模 2pi²，实际正体积是该模类在 (0,pi²) 的代表。
程序不通过小数选择这个代表；整数 level 的相位无需选择。
q1=4 sum sin(ell_j)sin(ell_opp)>0，Delta>0，故稳定根满足 |z|<1；全部 Z 同模，
所有 1-Z 具有正实部。这保证 Li₂/log 不触碰支割。Arg 的跳跃只改变 W 的整球体积。
正定相关矩阵域连通，驻点恒等式 exp(2z∂z L_e)=1 和正交例的 z∂z L_e=pi i
也排除了 q0=0；若 q0=0，则 z=0，左边的导数乘积应为零，矛盾。
这里正定性来自实际顶点，不能仅凭任意六个猜测边长的 det G>0 代替它。

## 3. 已完成的纯符号化简

一般约简器只使用有明确主支域的恒等式：

\[
\operatorname{Re}\operatorname{Li}_2(\bar z)=\operatorname{Re}\operatorname{Li}_2(z),\quad |z|<1;
\]
\[
\operatorname{Li}_2(z)+\operatorname{Li}_2(1-z)=\pi^2/6-\log z\log(1-z)
\quad (|z|<1,\ |1-z|<1);
\]
\[
\operatorname{Li}_2(z)+\operatorname{Li}_2(-z)=\tfrac12\operatorname{Li}_2(z^2),\quad |z|<1.
\]

未知域不强行约简；不无条件合并复对数，不用有理拟合。
这不是完备的 Li₂ 恒等式判定器；剩余项原样返回。

### 3.1 正交四面体

输入 (e1,e3,e1)，累计顶点为 (e0,e1,e2,e3)。所有 ell=pi/2，
q0=-4-4i、q1=12、q2=-4+4i、z=(1+i)/2。
四个正项的 Z 都为 z，四个负项都为共轭 z，故 Li₂ 实部完全抵消。
Sigma_j=-i pi、Arg(-q2)=-pi/4；于是

\[
W=-3\pi^2/8+\pi^2/4+3\pi^2/4-\pi^2/2=\pi^2/8,
\qquad \omega_k=e^{-ik\pi/8}.
\]

该结果既由代码符号约简获得，也由正交卦限占 S³ 的 1/16 独立确定。

### 3.2 双圆弧族的 Li₂ 微分恒等式

对于边长 (alpha,pi/2,pi/2,beta,pi/2,pi/2)，0<alpha,beta<pi，设
p=e^(i alpha)、q=e^(i beta)、S=pq+p+q-1、gamma=Arg S、eta=gamma-(alpha+beta)/2。
写 S=2e^(i(alpha+beta)/2)[cos((alpha-beta)/2)+i sin((alpha+beta)/2)]，
可得 0<eta<pi/2，并且 0<gamma<pi。此时 z=-2/S，八项归并为

\[
F=\operatorname{Li}_2(z)+\operatorname{Li}_2(-pqz)
-\operatorname{Li}_2(-pz)-\operatorname{Li}_2(-qz).
\]

恒等分解为

\[
1-z=\frac{(p+1)(q+1)}S,\quad
1+pqz=-\frac{(p-1)(q-1)}S,\quad
1+pz=\frac{(p+1)(q-1)}S,\quad
1+qz=\frac{(p-1)(q+1)}S.
\]

记对应四个主支对数为 LU、LV、LX、LY。逐项主辐角给出
LU+LV-LX-LY=-i pi，LX-LV=log cot(alpha/2)+i pi/2，
LY-LV=log cot(beta/2)+i pi/2。
利用 dLi₂(w)=-log(1-w)dlog w，得到
dF=i pi dlog z+i(LX-LV)dalpha+i(LY-LV)dbeta，
故 dRe F=pi dgamma-pi(dalpha+dbeta)/2=pi d eta；
以正交例 eta=pi/4、Re F=0 固定积分常数：

\[
\operatorname{Re}F=\pi\eta-\pi^2/4.
\]

六个 Im Sigma 都为 -pi，Arg(-q2)=gamma-pi，代回 W，所有额外项抵消，
得到 W=alpha beta/2。负 alpha 通过取向恢复 V_or=alpha beta/2。
独立几何参数化
x=(cos rho cos u,cos rho sin u,sin rho cos v,sin rho sin v)
给出体积形式 cos rho sin rho du∧d rho∧dv，同样积分为 alpha beta/2。
所以 (pi/3,pi/4) 与 (-pi/3,pi/4) 分别给出 ±pi²/24，六条边相同但相位互逆。

## 4. 2T：所有十二种非退化类型的实际体积

2T 的点为 {±e_j} 以及所有 h=(±1,±1,±1,±1)/2。
枚举全部 C(24,4)=10626 个无序集合，其中 5586 个秩不足 4，5040 个非退化。
对六个内积在 24 种顶点排列下取最小键，恰有 12 类。
可逆顶点矩阵具有相同 Gram 当且仅当由 O(4) 等距对应，所以这是无向几何分类，
不是群共轭类或仅凭六条边排序所得的伪分类。

下表 K=2(G01,G02,G03,G12,G13,G23) 的置换规范键。

| K | 集合重数 | 室数 m | Vol/pi² |
|---|---:|---:|---:|
| (-1,-1,-1,0,0,0) | 192 | 270 | 15/32 |
| (-1,-1,-1,0,1,1) | 576 | 84 | 7/48 |
| (-1,-1,0,0,-1,1) | 576 | 120 | 5/24 |
| (-1,-1,0,1,-1,-1) | 288 | 204 | 17/48 |
| (-1,-1,0,1,1,1) | 576 | 36 | 1/16 |
| (-1,-1,1,0,0,0) | 576 | 66 | 11/96 |
| (-1,-1,1,1,-1,0) | 576 | 60 | 5/48 |
| (-1,0,0,1,1,0) | 576 | 30 | 5/96 |
| (-1,0,1,1,0,1) | 576 | 24 | 1/24 |
| (0,0,0,0,0,0) | 48 | 72 | 1/8 |
| (0,0,1,0,1,1) | 192 | 18 | 1/32 |
| (0,1,1,1,1,1) | 288 | 12 | 1/48 |

### 4.1 为什么室计数就是实际体积

F4 根为短根 {±e_i}、所有半符号根，以及长根 {±e_i±e_j}。
它们给出 24 个不同中心超平面：4 个坐标面、12 个双坐标面、8 个全符号面。
带符号置换与矩阵

\[
H=\tfrac12\begin{pmatrix}1&1&1&1\\1&1&-1&-1\\1&-1&1&-1\\1&-1&-1&1\end{pmatrix}
\]

生成保持此根集的正交群 W。它在 24 个短根上传递；e0 稳定子必置换其余六个正交短根，
恰为 B3，阶 48。因此 |W|=24·48=1152。

任意三个独立 Hurwitz 顶点所张成的面都是这些超平面之一：先以 W 把第一个送到 e0；
其余为轴或半符号点。轴/轴给轴法向，轴/半符号点给双坐标法向；
两个半符号点在剩余三维投影为立方体顶点，其独立叉积只有两个等绝对值非零分量。
再经 W 变回，仍是 F4 根法向。

基本室可取 x1>x2>x3>0 且 x0>x1+x2+x3，其内部不被任何上述超平面穿过。
四个 inward 法向为 e1-e2、e2-e3、e3、(e0-e1-e2-e3)/2。
带根长颜色的角图为 3–4–3，颜色 long/long/short/short，唯一保持颜色与邻接的图自同构是恒等。
故 W 中室稳定子平凡；反射跨越相邻墙又使室作用传递。因此恰有 1152 个全等室，
每室体积 2pi²/1152=pi²/576。

每个非退化 Hurwitz 四面体的四面边界都在墙上，故它是完整小室的并，边界测度为零。
这一步使有限点判定成为严格体积计算，而不是抽样估计。

### 4.2 不含积分或小数的直接公式

令 P 为 (8,3,2,1)、(7,4,3,2)、(6,5,4,1) 的所有带符号排列，共 1152 个整点。
第一个严格在基本室内；H(8,3,2,1)=(7,4,3,2)，
H(8,3,2,-1)=(6,5,4,1)。三个各有 384 个点的不交轨道穷尽 W 轨道，故每室恰有一个点。

\[
m(V)=\#\{p\in P:V^{-1}p\text{ 的四分量严格为正}\},\qquad
V_{\rm or}=\operatorname{sgn}(\det V)\frac{m(V)\pi^2}{576}.
\]

全部判定通过整数余子式完成。这是**同一个球面四面体**的体积，没有 beta，没有代表转换。
必须先尝试该类型的 Li₂ 代入；无法消去的项再由上述几何恒等式识别其实际体积。
全类型的代表与逐类型尝试记录由端到端 JSON 给出。

## 5. 2O 与 2I：完整类型参数和不可有理化的严格反例

2O=2T 加上全部 (±e_i±e_j)/sqrt(2)。
2I=2T 加上 (0,±1/2,±phi/2,±phi^(-1)/2) 的全部偶坐标置换，phi=(1+sqrt(5))/2。
坐标分别落在 Q(sqrt(2))、Q(sqrt(5))，程序以整数对表示，所有秩判定精确。

利用左乘等距性，把一个顶点移到 e0，然后枚举其余三个不同点；
Gram 仍在 24 种顶点排列下规范化。每个原无序集合有四种被选中的基点，
因此全部集合重数=锚定重数·|G|/4；原有序三输入重数=锚定重数·6。

| 群 | 类型数 | 锚定非退化集合 | 锚定退化集合 | 全部非退化四点集合 | 非退化三输入 |
|---|---:|---:|---:|---:|---:|
| 2T | 12 | 840 | 931 | 5040 | 5040 |
| 2O | 84 | 10080 | 6135 | 120960 | 60480 |
| 2I | 563 | 212520 | 61299 | 6375600 | 1275120 |

所有 84/563 行均给出代表、完整 Gram 键、重数、六个精确余弦和 det G；
共用第 2 节完整显式公式即可求相位，无需数值体积数据库。
JSON 的 `formula_template` 是顺序赋值的精确解析程序，每行 `formula_inputs` 给出全部代入参数。
PDF 内另外印出全部 84/563 个六位 Gram 参数键、锚定重数和解码字母表；
这些参数直接决定六个 c_j 和 det G，因此完整参数表不只存在于程序接口。
可通过 `--group 2I --type-index 0` 单独导出任意类型的展开式。
符号化简尝试使用明确的 CAS 时间预算；超时表示此次尝试未结束，不证明任何不可化简定理。

### 5.1 2O 的反例

取累计顶点 e0、e1、(e1+e2)/sqrt(2)、(e2+e3)/sqrt(2)，D=1/2。
对应输入为
g1=e1，g2=(e0+e3)/sqrt(2)，g3=(e0+e1-e2+e3)/2，均属于 2O。
后三点构成 S² 三角形，其角为 pi/4、pi-atan sqrt(2)、pi/2-atan sqrt(2)。
Girard 公式给面积 Omega=3pi/4-2atan sqrt(2)。
令 theta=Omega/2，则 tan theta=1/(3+sqrt(2)) 且 0<theta<pi/8。
与正交点 e0 作球面 join，体积形式为 sin²u du dArea，u∈[0,pi/2]，因此

\[
V_{2O}=\frac\pi2\arctan\frac1{3+\sqrt2},\qquad
\omega_k=\exp\left(-\frac{ik}{2}\arctan\frac1{3+\sqrt2}\right).
\]

设 t=tan theta，则 e^(2i theta)=(1+it)/(1-it)∈Q(sqrt(2),i)=Q(zeta8)。
该域的单位根仅为 mu8，而 0<2theta<pi/4 排除其所有元素。
故 theta/pi、V/pi² 都不是有理数，非零整数 k 的该相位也不是单位根。

### 5.2 2I 的反例

取累计顶点 e0、e1、e2、(0,1/2,phi/2,phi^(-1)/2)，D=phi^(-1)/2>0。
对应输入 e1、e3、(phi/2,phi^(-1)/2,0,-1/2) 均属于 2I。
空间三角形的半球面积正切为
det(A,B,C)/(1+A·B+B·C+C·A)=1/(4phi+1)，可由 Girard 公式及半角恒等式得到。
同一个球面 join 积分给出

\[
V_{2I}=\frac\pi2\arctan\frac1{4\varphi+1},\qquad
\omega_k=\exp\left(-\frac{ik}{2}\arctan\frac1{4\varphi+1}\right).
\]

e^(2i theta) 属于双二次域 Q(sqrt(5),i)。其单位根只有 mu4：
可能的四次或更低圆分域中，mu3/mu6 要求 sqrt(-3)，mu8 要求 sqrt(2)，mu12 要求 sqrt(3)；
mu5/mu10 给出循环四次域，不能等于此双二次域。更高阶的圆分次数也不能整除 4。
因 0<2theta<pi/2，它不是单位根，从而同样排除有理体积和非零 level 的单位根相位。
2O 的 mu8 判断也由同一低次数圆分域检查得到。

因此“2O/2I 通过坐标符号规律全部给出有理体积”不是可实现的正确目标。
这个结论是解析证明，不来自数值残差或拟合失败。
它不与有限群反常类具有有限阶矛盾：类的 N 次幂平凡只意味着代表的 N 次幂是某个余边界，
并不意味着原几何代表的每个逐点值都满足 omega^N=1。本轮计算没有借此更换代表。

## 6. 退化输入：原 E/F/T 的精确解析覆盖

全局接口先作相邻重复归一化，再精确判定零是否属于顶点凸包。
短边不安全时通过固定 J=e1 分段；不安全的面/三链按原规则锥接边界，
锥点取首个不落在相关顶点张成空间内的 a(1,t,t²,t³)，0≤t≤3M。
每个真子空间至多排除三个候选，保证有限终止。每面至多 6 项，每三链至多 24 项。

秩不足 4 的安全项积分为零，但保留其链级记录；满秩项按正范数单位化后用第 2 节公式。
最终是每项原系数、定向符号与 W 的有限和，再取 exp(-ik sum/pi)。
全局链边界相容性与原构造一致，因此五边形由闭三链整球体积证明，
不是由这次枚举若干四元组声称成立。

关键限制：锥点通常不在 E 群内，所以不能把它们强行映射到 12/84/563 类型表。
2T 十二个有理体积也不等于已经证明了所有退化 E/F/T 单形都有有理体积。
`expand=False` 仅查看精确链，有未展开非零体积项时 `phase=None`，不会冒充已求得相位。
默认展开返回完整主支特殊函数表达式；`reduce=True` 可再请求较昂贵的 CAS 化简。

## 7. 入口与可复跑证据

```python
from su2_symbolic import group_formula
from sympy import pi
r = group_formula((0,1,0,0),(0,0,0,1),(0,1,0,0))
assert r['raw_volume'] == pi**2/8

from su2_exact_ade import global_formula, classify_group
exact_phase = global_formula('2I',(-1,0,0,0),(-1,0,0,0),(-1,0,0,0))
types = classify_group('2I')
```

```bash
.venv/bin/python -m unittest tests.test_su2_symbolic tests.test_2t_geometry tests.test_exact_ade tests.test_exact_certificates tests.test_symbolic_cli -v
.venv/bin/python -m scripts.check_symbolic_2t --symbolic-trials all --symbolic-budget-seconds 4 --jobs 4 --output results/SU2_symbolic_ADE_exact.json
.venv/bin/python -m scripts.check_symbolic_2t --group 2O --type-index 0 --output results/example_2O_type0.json
.venv/bin/python scripts/build_complete_report.py --engine tectonic
```

新模块拒绝浮点群坐标、非单位元、非整数 level 和非法群成员。
历史含 mpmath 的测试只作为旧接口兼容性回归，不参与本稿的符号化简、类型赋值或无理性证明。
新数学模块没有数值求值入口。新结论是可提交专家审阅的研究结果，不声称已通过外部同行评议。

最终回归为 70 项通过。659 类型的四秒符号尝试结果（已消去 Li₂ / 明确剩余项 / 预算内未完成）
分别为：2T 的 3/9/0，2O 的 17/42/25，2I 的 33/40/490。
这里“消去 Li₂”并不保证 log/Arg 也已消去，更不等价于有理体积。
完整 PDF 45 页，印有所有 84/563 个参数键；详见 [复核和验收记录](REVIEW_AND_VALIDATION.md)。
