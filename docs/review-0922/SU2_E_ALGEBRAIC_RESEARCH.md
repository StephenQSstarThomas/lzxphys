# E 型后续研究：精确轨道计数公式与代表转换

日期：2026-09-22。性质：本轮新推导及其实现、经独立数学与源码复核的专家审阅稿；代理辅助复核不等于外部同行评审。

## 0. 回应本轮 review 的范围

是，本轮人类 review 要求的不只是“E 型数值通过”，而是从三个群元素出发、能独立读懂的相位公式及其依据。必要后续项包括：公式的全部分支、闭性、归一化、反常类的识别、退化输入，以及可核查的验证记录。百万行相位表、全部四元组穷举、类似 D 型八行表的最短表达，则不应擅自当成这些要求的同义词。

本补充的实质进展是一个**不含体积积分、反三角函数或二重对数的精确代数代表**。它保留既有的有限 E/F/T 链，但将其体积替换为有向轨道计数。以下给出公式、证明、与原几何代表的显式余边界关系，以及逐个可检查的例子。

重要限制：新代表与原几何代表一般不逐点相等。用户已明确批准“保留原公式，并证明转换关系”的方案。因此保留 `phase_from_exact_quaternions` 不变，另在 `scripts/ade_count.py` 提供 `count_details`、`algebraic_phase` 和 `gauge_details`，不静默替换旧代表。

## 1. 输入、字段和有限链

写四元数为 $q=(a,\mathbf u)$，对应 $aI+i\mathbf u\cdot\boldsymbol\sigma$。乘法固定为

\[
(a,\mathbf u)(b,\mathbf v)
=(ab-\mathbf u\cdot\mathbf v,\ a\mathbf v+b\mathbf u-\mathbf u\times\mathbf v).
\]

令 $e_0,e_1,e_2,e_3$ 为标准坐标基，$\varphi=(1+\sqrt5)/2$。三个群的**全部**元素为

\[
\begin{aligned}
\mathcal Q_6={}&\{\pm e_j:0\le j\le3\}
\cup\{(\epsilon_0,\epsilon_1,\epsilon_2,\epsilon_3)/2:\epsilon_j=\pm1\},\\
\mathcal Q_7={}&\mathcal Q_6\cup
\{(\epsilon e_i+\epsilon'e_j)/\sqrt2:i<j,\ \epsilon,\epsilon'=\pm1\},\\
\mathcal Q_8={}&\mathcal Q_6\cup
\{\rho(0,\epsilon_1/2,\epsilon_2\varphi/2,\epsilon_3\varphi^{-1}/2):
\rho\in A_4,\ \epsilon_j=\pm1\}.
\end{aligned}
\]

这里 $A_4$ 表示偶坐标置换，$\Gamma=\mathcal Q_m$，群阶 $N=24,48,120$，对应字段为 $\mathbb Q,\mathbb Q(\sqrt2),\mathbb Q(\sqrt5)$。

给定 $g_1,g_2,g_3\in\Gamma$，令

\[
(p_0,p_1,p_2,p_3)=(e_0,g_1,g_1g_2,g_1g_2g_3).
\]

使用 [上一轮完整公式的 §3.2](SU2_HUMAN_REVIEW_RESPONSE.md) 已定义的 E/F/T 链，展开

\[
T(p_0,p_1,p_2,p_3)=\sum_{\nu=1}^{M\le24}c_\nu
[v_{\nu0},v_{\nu1},v_{\nu2},v_{\nu3}].
\tag{1}
\]

为避免把这个记号变成黑箱，这里复述全部选择规则。记 $S(v_0,\dots,v_j)$ 为不存在非零全非负系数 $\lambda_i$ 使 $\sum_i\lambda_i v_i=0$。相邻重复顶点给零链；固定 $J=e_1$，其余分支为

\[
E(a,b)=\begin{cases}[a,b],&S(a,b),\\[a,aJ]+[aJ,b],&\text{否则};\end{cases}
\]

对含 $r$ 项的链 $Z=\sum_\mu d_\mu[w_{\mu0},\dots,w_{\mu j}]$，取首个 $t\in\{0,\dots,3r\}$ 使 $u=a(1,t,t^2,t^3)$ 不在任何一项顶点的线性张成中，定义

\[
C_a Z=\sum_\mu d_\mu[u,w_{\mu0},\dots,w_{\mu j}],\qquad C_a0=0.
\]

每个待避真子空间至多排除三个候选，故搜索必终止。继而

\[
F(a,b,c)=\begin{cases}[a,b,c],&S(a,b,c),\\
C_a\bigl(E(b,c)-E(a,c)+E(a,b)\bigr),&\text{否则},\end{cases}
\]

\[
T(a,b,c,d)=\begin{cases}[a,b,c,d],&S(a,b,c,d),\\
C_a\bigl(F(b,c,d)-F(a,c,d)+F(a,b,d)-F(a,b,c)\bigr),&\text{否则}.\end{cases}
\]

面链至多 6 项、19 个候选；三维链至多 24 项、73 个候选。所有输出单形安全。降秩项保留在链级，只在三维计数及积分时贡献零。单形是单位化顶点的径向参数化；存储时可使用任意正比例射线，不能把负比例视作相同顶点，也不能将倒序单形直接当作负单形约去。

## 2. 新公式：四阶行列式和精确正负号

对 $D_\nu=\det(v_{\nu0},\dots,v_{\nu3})\ne0$，$h\in\Gamma$，$j,s\in\{0,1,2,3\}$，定义

\[
d_{\nu hjs}=
\det(v_{\nu0},\dots,v_{\nu,j-1},he_s,v_{\nu,j+1},\dots,v_{\nu3}).
\tag{2}
\]

令 $\operatorname{lsign}(a_0,a_1,a_2,a_3)$ 表示从左向右第一个非零元素的符号。本处每组系数必不全为零：对应的非零 Cramer 线性泛函不可能同时湮灭基 $he_0,\dots,he_3$。令

\[
I_{\nu h}=\prod_{j=0}^{3}
\mathbf1\!\left\{\operatorname{lsign}(d_{\nu hj0},d_{\nu hj1},d_{\nu hj2},d_{\nu hj3})
=\operatorname{sgn}D_\nu\right\}.
\tag{3}
\]

于是整数和及相位为

\[
\boxed{A_\Gamma(g_1,g_2,g_3)=
\sum_{\nu:D_\nu\ne0}c_\nu\operatorname{sgn}D_\nu
\sum_{h\in\Gamma}I_{\nu h}\in\mathbb Z,}
\tag{4}
\]

\[
\boxed{\widehat\omega_{\Gamma,k}(g_1,g_2,g_3)
=\exp\!\left(-\frac{2\pi i k}{N}A_\Gamma(g_1,g_2,g_3)\right),
\qquad k\in\mathbb Z.}
\tag{5}
\]

三个 E 型公式分别取 $(\Gamma,N)=(\mathcal Q_6,24),(\mathcal Q_7,48),(\mathcal Q_8,120)$。式 (1)–(5) 的所有输入和有限选择均已给出。每次至多 $24N$ 个“单形—轨道点”判定；这是 576、1152、2880 个判定的统一上界，不是同样次数的单次算术运算。

等价且便于实现的形式：$V_\nu$ 以顶点为列，$L_h$ 是左乘 $h$ 的四阶矩阵。$I_{\nu h}=1$ 当且仅当 $V_\nu^{-1}L_h$ 的四行各自第一个非零元素均为正。

所有算术留在原字段。对 $a+b\sqrt d$ 的符号，当 $a,b$ 同号或其中之一为零时直接判断；异号时比较 $a^2$ 与 $db^2$。因 $d=2,5$ 非平方，非零有理 $a,b$ 不会产生相等歧义。无需浮点阈值、选择一个“足够小”的机器数，或猜测相位的有理倍数。

这是一条精确、有界的代数求和公式，**不是已经导出的 E 型专属短分支表**。

## 3. 全部三元组的闭性证明

### 3.1 无穷小只是确定性的符号规则

考察 $r(\varepsilon)=(1,\varepsilon,\varepsilon^2,\varepsilon^3)$，$\varepsilon>0$。Cramer 法则给出

\[
\bigl(V_\nu^{-1}L_h r(\varepsilon)\bigr)_j
=D_\nu^{-1}\sum_{s=0}^{3}d_{\nu hjs}\varepsilon^s.
\]

充分小的正 $\varepsilon$ 下，其符号正是式 (3)。群有限、每条链有限，因此全部 $\Gamma^4$ 输入产生的链、面、降秩项及其群平移只涉及有限个真线性子空间。每个真子空间可以用一个非零线性泛函排除；该泛函沿 $r(\varepsilon)$ 是非零、次数至多三的多项式。所以存在一个共同的充分小正实数，使全部符号规则同时成立，且探测点避开所有低维像。

这给出了普通实球面中的解释，而不是使用未经说明的非阿基米德拓扑。公式求值本身仍只取首个非零系数的符号，不需要找出这个实数。

### 3.2 轨道计数是局部次数

对安全满秩径向单形 $Q=[v_0,v_1,v_2,v_3]$，单位点 $y$ 在其内部，当且仅当 $V^{-1}y$ 四个坐标全正。该单形在内部的有向局部重数为 $\operatorname{sgn}\det V$。因此式 (4) 是链 $T$ 在 $N$ 个探测点 $h\,r(\varepsilon)/\|r(\varepsilon)\|$ 上的有向重数之和。

左乘保向，且置换这一轨道，所以 $A_\Gamma$ 的齐次版本左乘不变。这里必须平均整个群轨道；只数一个固定点通常不具有所需等变性。

对任意闭三链 $Z$ 和避开各项边界及低维像的 $y$，局部次数定理给出

\[
\chi_y(Z)=\deg Z=\frac1{2\pi^2}\int_ZdV\in\mathbb Z.
\tag{6}
\]

E/F/T 使用规范化奇异链。若将规范化闭链提升为普通奇异闭链，标准规范化分裂仅加入面的退化项；它们的像仍在上述低维集合中，对三维积分及所选点的局部重数均为零。因此式 (6) 同样适用。

### 3.3 五边形与归一化

边界恒等式 $\partial T=\delta F$ 给出 $\partial(\delta T)=0$。于是

\[
\delta A_\Gamma=N\deg(\delta T)\in N\mathbb Z.
\tag{7}
\]

非齐次写法即

\[
A(b,c,d)-A(ab,c,d)+A(a,bc,d)-A(a,b,cd)+A(a,b,c)
\equiv0\pmod N.
\]

式 (5) 因而满足完整五边形。任一输入为单位元时，累计顶点出现相邻重复，$T=0$、$A=0$，相位严格等于 1。这一证明覆盖退化和反足输入，并非由若干抽样数值推断。

## 4. 与原几何相位的显式转换

这里必须区分“同一个反常类”和“同一个逐点函数”。令

\[
t(g,h,l)=\frac{1}{2\pi^2}\int_{T(1,g,gh,ghl)}dV,
\qquad \tau(g,h,l)=\frac{A_\Gamma(g,h,l)}N.
\]

既有边长公式给出 $t$ 的实际有限求值。由式 (6)–(7)，$f=\tau-t$ 是**实数值**三余循环：$\delta f=0$，而不只是模整数为零。定义明确的二余链

\[
B(g,h)=\frac1N\sum_{x\in\Gamma}f(x,g,h),
\qquad \beta_k(g,h)=\exp[-2\pi i k B(g,h)].
\tag{8}
\]

对 $\delta f(x,g,h,l)=0$ 按 $x$ 求和，并用 $x\mapsto xg$ 置换群，得到

\[
f(g,h,l)=B(h,l)-B(gh,l)+B(g,hl)-B(g,h)=\delta B(g,h,l).
\]

因此有**逐点成立的代表转换恒等式**

\[
\boxed{\widehat\omega_{\Gamma,k}
=\omega^{\rm geom}_{A,k}\,\delta\beta_k,\qquad
\delta\beta_k(g,h,l)=
\frac{\beta_k(h,l)\beta_k(g,hl)}{\beta_k(gh,l)\beta_k(g,h)}.}
\tag{9}
\]

$B$ 和 $\beta$ 均归一化。式 (8) 是实际有限和，不只是声称某个二余链存在。计算新相位无需先计算 $B$；仅在要比较两个代表时才需原几何体积。

## 5. 类阶：不以小数或循环子群抽查代替

还可以把类识别的拓扑步骤明确写出。$\Gamma$ 左乘自由且保向地作用于单位四元数球面。相应有向实四维表示的球面丛为

\[
S^3\longrightarrow E\Gamma\times_\Gamma S^3\longrightarrow B\Gamma,
\qquad E\Gamma\times_\Gamma S^3\simeq S^3/\Gamma.
\]

商空间是定向闭三流形，故 Gysin 序列包含

\[
H^3(S^3/\Gamma;\mathbb Z)\cong\mathbb Z
\xrightarrow{\ \times N\ }H^0(B\Gamma;\mathbb Z)\cong\mathbb Z
\xrightarrow{\ \smile e\ }H^4(B\Gamma;\mathbb Z)\longrightarrow0.
\tag{10}
\]

第一箭头是沿纤维积分；纤维包含映射在上述同伦等价下是覆盖 $S^3\to S^3/\Gamma$，其次数为 $N$。于是 $H^4(B\Gamma;\mathbb Z)=\mathbb Z/N$，由 Euler 类 $e$ 生成。

E/F/T 依次填充 bar 单形的零至三维边界；四维单形边界的球面次数 $\deg(\delta T)$ 正是这一球面丛的首个障碍余循环，即定向选择下的 $e$。链级解释使用 $S^3$ 的二连通性及 $\pi_3(S^3)\cong H_3(S^3)$；不要求将每个有限链误认为一个单独的非退化四面体。

通过 $0\to\mathbb Z\to\mathbb R\to U(1)\to0$ 的连接同态，$\exp(-2\pi i k t)$ 映到 $-k[\deg(\delta T)]$。有限群的正次实上同调为零，故连接同态是同构。配合式 (9)，两代表的类阶均为

\[
\boxed{\operatorname{ord}[\widehat\omega_{\Gamma,k}]
=\operatorname{ord}[\omega^{\rm geom}_{A,k}]
=\frac{N}{\gcd(k,N)}.}
\tag{11}
\]

这也与 [Epa–Ganter，定理 1.1、命题 4.1](https://arxiv.org/pdf/1605.09192v1) 对有限球面子群的结论一致。整体正负号不影响阶；与源 A 的正负约定仍沿用仓库已有定义，不能声称重新检查了当前缺失的原附件。

注意：$C_4$ 配对仅能校准该限制。三类群的循环子群阶的最小公倍数分别只有 12、24、60；单靠这些循环限制不能识别全部 24、48、120 阶，尤其不能遗漏二主部分的额外因子 2。

## 6. 三个明确输入，逐个列出贡献者

共同取 $g_1=e_1,g_2=e_3$，沿用原文的三个 $g_3$。每个例子的链都是一个正取向直接四面体，故整数 $A$ 就是贡献轨道点数。

| 群 | $g_3$ | 使 $I_{\nu h}=1$ 的全部 $h$ | $A$ | 新代数相位 |
|---|---|---|---:|---|
| $2T$ | $(1,1,-1,-1)/2$ | $e_0$ | 1 | $e^{-i\pi k/12}$ |
| $2O$ | $(1,1,0,0)/\sqrt2$ | $e_0,(1,0,1,0)/\sqrt2$ | 2 | $e^{-i\pi k/12}$ |
| $2I$ | $(\varphi/2,\varphi^{-1}/2,0,-1/2)$ | $e_0,(1,\varphi^{-1},\varphi,0)/2,(\varphi,1,\varphi^{-1},0)/2,(\varphi^{-1},\varphi,1,0)/2$ | 4 | $e^{-i\pi k/15}$ |

此表在本轮由主审和两路独立计算分别复现，采用精确 SymPy 代数符号判定，不使用浮点阈值。原几何代表对应的相位仍分别是

\[
e^{-i\pi k/32},\qquad e^{-i\pi k/16},\qquad
\exp\!\left[-\frac{ik}{2}\arctan\frac1{4\varphi+1}\right].
\]

两行结果不同是式 (9) 的必要体现，不是数值误差。不能把新表填回原几何接口的预期值。

### 为什么 E8 原例甚至不是单位根

设 $u=(4\varphi+1)^{-1}$，$\theta=\arctan u$。则

\[
e^{2i\theta}=\frac{1+iu}{1-iu}\in K=\mathbb Q(\sqrt5,i),\qquad 0<\theta<\pi/4.
\]

$K$ 的单位根只有 $\{\pm1,\pm i\}$。证明：若含原始 $m$ 次根，则 $\phi(m)$ 整除 4，故只需检查 $m=1,2,3,4,5,6,8,10,12$。$K$ 的三个二次子域恰为 $\mathbb Q(\sqrt5),\mathbb Q(i),\mathbb Q(\sqrt{-5})$，排除含 $\sqrt{-3},\sqrt2,\sqrt3$ 的情形；$\mathbb Q(\zeta_5)$ 是循环四次扩张，也不是 Galois 群为 $C_2\times C_2$ 的 $K$。因此只剩 $m=1,2,4$。

上述开区间又排除 $e^{2i\theta}$ 等于这四个根，故 $\theta/\pi$ 无理，任何非零整数 $k$ 的原几何相位 $e^{-ik\theta/2}$ 均不是单位根。

这严格排除了“原几何相位逐点等于 $\mu_{120}$ 值表”的方案；并不排除含反三角函数的其他短表达。

## 7. 本轮实际验证与交付状态

本轮实际执行并已持久化的检查：

- 三群全部 24、48、120 个元素：无重复、精确范数为 1、精确共轭逆元属于群且乘积为单位元。另重跑 `check_ade_review`，确认三个列表的全部两两乘积闭包；该脚本同时重跑 Q8 的 4096 个几何五边形，最大残差约 $1.13746\times10^{-86}$。
- 每群 6 个累计顶点样例，共 18 个：逐项验证 $\partial F=\delta E$、$\partial T=\delta F$、输出单形安全、左乘最后一个列表元素后的精确链等变性。前五类为 $(1,e_1,e_2,p_3)$、$(1,x,1,x)$、$(1,-1,x,1)$、$(1,-1,1,-1)$、$(1,1,x,1)$，其中 $x$ 为上表该群的 $g_3$，$p_3=e_2x$；链项数分别为 1、1、13、8、0。第六类为 $(1,r,r^2,x)$，$r=(-1,1,1,1)/2$，其中 $1+r+r^2=0$ 且没有反足点。
- 每群各一个直接分支和锥分支的正射线缩放检查；二次字段平方关系及黄金分割关系检查。
- 式 (2)–(4) 对上表三个精确输入的计数分别为 1、2、4；所有贡献者完整列在上表，不仅报告“通过”。
- 单元套件共 39 项，另覆盖中心相位、三个位置的单位元归一化、非成员/浮点/未知字段输入拒绝、正根嵌入、负 level 与 $10^{120}+1$ 的精确模运算，以及 β 转换的符号。
- 每个 E 群 16 组确定性四元组，共 48 个整数模 $N$ 五边形：全部余数严格为 0。每群 $C_4$ bar 循环的指数分别为 6/24、12/48、30/120，均严格给出 $+i$。
- 每群上表一个三元组的完整 β 转换检查：计算四个二余链值，每个包含全部 $N$ 个平均项，共 $4(24+48+120)=768$ 个平均项。45 位输出下三群残差依次约为 $4.54848\times10^{-46}$、$3.15172\times10^{-46}$、$6.11528\times10^{-46}$，验收阈值为 $10^{-39}$。

完整记录为 [SU2_E_algebraic_validation.json](SU2_E_algebraic_validation.json)：包括每个五边形的全部五个整数计数、模余数、精确输入、所有单形与贡献轨道元素，以及 β 的全部平均项。它不是 E 型四元组穷举；浮点转换检查也不是区间误差证书。

复现（从仓库根目录、已安装 `.[validation]` 的环境运行）：

```bash
python -m unittest discover -s tests -v
python -m scripts.check_e_algebraic
python -m scripts.check_ade_review
python scripts/build_complete_report.py
```

`check_e_algebraic --skip-gauge` 只做快速精确计数检查，报告会明确标记未运行 β 检查，不能把它作为完整运行的替代证据。

独立源码审计另外对 36 个有理/二次字段矩阵比较了实现与 SymPy 的伴随矩阵、对 180 个二次字段值比较了符号，并检查负叉积左乘矩阵与缓存隔离。主审修复了一个审计中发现的输入缺陷：有理字段曾把未知符号 `s0` 当作零；现明确拒绝，回归测试先失败后通过。

本轮交付可进入专家审阅：完整公式、全部分支定义、闭性与归一化证明、显式代表转换、类阶证明、精确实现及可复现验证。**仍未交付**的是类似 D 型八行表的 E 型专属最短表，以及 E 型所有四元组的穷举；这两项不应与当前有界精确公式混为一谈。任意 $n$ 的 D 表证明属于另一项遗留任务，不由本轮 E 型结果代为解决。

另一条可能路线是有限群周期分解到 bar 分解的显式比较映射。[Tomoda–Zvengrowski，§4](https://arxiv.org/pdf/0904.1876) 给出了相关周期分解；但把它变成 E8 的易读群元素三余循环仍需进一步推导。本轮优先使用上面的有界计数路线，不宣称已完成该替代项目。
