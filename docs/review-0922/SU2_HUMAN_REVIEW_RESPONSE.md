# 人类 review 修正：三个群元素的边长公式与 E 型显式相位

本轮以用户提出的两个缺口为验收标准。当前 checkout 中没有名为 09-21 的目录或原始任务附件；可核查的推进材料是 `docs/review-0922/`、两份 LaTeX 文稿和对应源码。此前“已完成”的交接不能代替下面的实际公式。

核查结论：边长路线已经在复核附录中出现，但 ADE 主答案仍重列二面角公式，没有把三个群元素的六个迹连到最终答案。E 型已有一个几何求值实现，但主文档只指向函数，旧 JSON 每群只留两组检查的最大残差和前三项相位；这不构成可独立阅读、逐项复算的数学交付。

## 1. 三个群元素直接给出六条边长

固定单位半径的 $S^3\simeq SU(2)$，写
\[
g=q_0I+i(q_1\sigma_1+q_2\sigma_2+q_3\sigma_3),\qquad q\in S^3.
\]
本坐标下的乘法是
\[
(a,\mathbf u)(b,\mathbf v)
=(ab-\mathbf u\cdot\mathbf v,\ a\mathbf v+b\mathbf u-\mathbf u\times\mathbf v).
\tag{1}
\]
取累计顶点
\[
p_0=1,\quad p_1=g_1,\quad p_2=g_1g_2,\quad p_3=g_1g_2g_3.
\]
因为 $q(p)\cdot q(r)=\frac12\operatorname{Tr}(p^{-1}r)$，短测地线长度为
\[
\ell_{ij}=\arccos\!\left(\tfrac12\operatorname{Tr}(p_i^{-1}p_j)\right)\in[0,\pi].
\]
因此，按后续公式固定的编号，**输入只需三个群元素**：
\[
\boxed{
\begin{aligned}
\ell_1=\ell_{23}&=\arccos\tfrac12\operatorname{Tr}g_3,&
\ell_2=\ell_{13}&=\arccos\tfrac12\operatorname{Tr}(g_2g_3),\\
\ell_3=\ell_{03}&=\arccos\tfrac12\operatorname{Tr}(g_1g_2g_3),&
\ell_4=\ell_{01}&=\arccos\tfrac12\operatorname{Tr}g_1,\\
\ell_5=\ell_{02}&=\arccos\tfrac12\operatorname{Tr}(g_1g_2),&
\ell_6=\ell_{12}&=\arccos\tfrac12\operatorname{Tr}g_2.
\end{aligned}}
\tag{2}
\]
这里是球面弧长，不是 $\mathbb R^4$ 中的弦长。Cayley–Menger 公式计算欧氏四面体体积，不能替代球面反常体积。

## 2. 只用六边长的完整体积表达

下面是 [Murakami 定理 1.2](https://arxiv.org/pdf/1011.2584v4) 的边长式，已展开其中所有固定 $z$ 偏导。这里不求逆 Gram 矩阵，也不先算原四面体的二面角。

令 $\bar j=j+3\pmod6\in\{1,\ldots,6\}$，定义
\[
\tau_j=\pi-\ell_{\bar j},\qquad a_j=-e^{-i\ell_{\bar j}}=e^{i\tau_j},
\]
\[
\Gamma=\begin{pmatrix}
1&\cos\ell_4&\cos\ell_5&\cos\ell_3\\
\cos\ell_4&1&\cos\ell_6&\cos\ell_2\\
\cos\ell_5&\cos\ell_6&1&\cos\ell_1\\
\cos\ell_3&\cos\ell_2&\cos\ell_1&1
\end{pmatrix}.
\tag{3}
\]
对于非退化累计四面体，$\det\Gamma>0$。所有辅助量均由边长给出：
\[
\begin{aligned}
r_0={}&a_1a_4+a_2a_5+a_3a_6+a_1a_2a_6+a_1a_3a_5
+a_2a_3a_4+a_4a_5a_6+\prod_{j=1}^6a_j,\\
r_1={}&4\sum_{j=1}^3\sin\tau_j\sin\tau_{j+3},\qquad r_2=\overline{r_0},\\
z={}&\frac{-2r_0}{r_1+4\sqrt{\det\Gamma}},\qquad |z|<1.
\end{aligned}
\tag{4}
\]
此处的判别式为 $r_1^2-4|r_0|^2=16\det\Gamma$，其中是**顶点** Gram；不能沿用原二面角路线的逆 Gram 行列式。

记
\[
\mathcal P=(1245,1346,2356),\quad
\mathcal N=(123,156,246,345),\quad A_S=\prod_{r\in S}a_r.
\]
一个编号串表示它包含的指标集合。八项二重对数明确为
\[
L_e=\frac12\left[\operatorname{Li}_2(z)
+\sum_{S\in\mathcal P}\operatorname{Li}_2(z/A_S)
-\sum_{S\in\mathcal N}\operatorname{Li}_2(-z/A_S)
-\sum_{j=1}^3\tau_j\tau_{j+3}\right].
\tag{5}
\]
令 $P_r=\Im\log(1-z/A_{\mathcal P_r})$、$N_r=\Im\log(1+z/A_{\mathcal N_r})$。六个实数是
\[
\begin{aligned}
d_1&=(P_1+P_2-N_3-N_4+\tau_1)/2,\\
d_2&=(P_1+P_3-N_2-N_4+\tau_2)/2,\\
d_3&=(P_2+P_3-N_2-N_3+\tau_3)/2,\\
d_4&=(P_1+P_2-N_1-N_2+\tau_4)/2,\\
d_5&=(P_1+P_3-N_1-N_3+\tau_5)/2,\\
d_6&=(P_2+P_3-N_1-N_4+\tau_6)/2.
\end{aligned}
\tag{6}
\]
定义完全展开的边长函数
\[
\boxed{\mathcal B(\ell_1,\ldots,\ell_6)
=\left[\Re L_e-\pi\operatorname{Arg}(-r_2)
-\sum_{j=1}^6\ell_jd_j-\frac{\pi^2}{2}\right]_{2\pi^2}.}
\tag{7}
\]

所有 $\operatorname{Li}_2$、$\log$、Arg 取主支，平方根取正实根；$[x]_{2\pi^2}=x-2\pi^2\lfloor x/(2\pi^2)\rfloor$。非退化凸球面四面体的结果属于 $(0,\pi^2)$。式 (6) 的来源是固定 $z$ 时 $\partial_{\ell_j}\log a_{\bar j}=-i$，配合 $d\operatorname{Li}_2(x)=-\log(1-x)d\log x$；因此没有留下未计算的导数。

六条边长对镜像不变，故还必须从群元保留
\[
D=\det(q(p_0),q(p_1),q(p_2),q(p_3)).
\]
最终一般位置答案为
\[
\boxed{V(g_1,g_2,g_3)=\operatorname{sgn}(D)\mathcal B(\ell(g_1,g_2,g_3)),\qquad
\omega_{A,k}=\exp[-ikV/\pi],\quad D\ne0.}
\tag{8}
\]
这已经把三个群元素与体积直接连接。$D=0$ 不意味着全局相位等于 1；下一节给出退化输入的有限公式。

## 3. E6、E7、E8 的群元素相位公式

### 3.1 输入的精确集合

令 $e_0,e_1,e_2,e_3$ 为四维坐标基，$\varphi=(1+\sqrt5)/2$。三个集合为
\[
\begin{aligned}
\mathcal Q_6={}&\{\pm e_r:0\le r\le3\}
\ \cup\ \{(\epsilon_0,\epsilon_1,\epsilon_2,\epsilon_3)/2:\epsilon_r=\pm1\},\\
\mathcal Q_7={}&\mathcal Q_6\ \cup\
\{(\epsilon e_r+\epsilon'e_s)/\sqrt2:r<s,\ \epsilon,\epsilon'=\pm1\},\\
\mathcal Q_8={}&\mathcal Q_6\ \cup\
\{\rho(0,\epsilon_1/2,\epsilon_2\varphi/2,\epsilon_3\varphi^{-1}/2):
\rho\in A_4,\ \epsilon_r=\pm1\}.
\end{aligned}
\tag{9}
\]

它们分别是 $2T,2O,2I$，阶为 24、48、120。任取 $q_1,q_2,q_3\in\mathcal Q_m$，用式 (1) 乘法构造 $p=(e_0,q_1,q_1q_2,q_1q_2q_3)$。计算在 $\mathbb Q$、$\mathbb Q(\sqrt2)$、$\mathbb Q(\sqrt5)$ 内完成；根取正根。

### 3.2 覆盖所有三元组的有限系数规则

以下给出式 (10) 所有求和项的定义，读者不需要运行 Python 才能知道相位是什么。符号 $[v_0,\ldots,v_j]$ 表示单位化顶点的径向有向单形；实际计算可保留未归一化的正射线。令
\[
S(v_0,\ldots,v_j)\iff 0\notin\operatorname{conv}\{v_0,\ldots,v_j\}.
\]
它等价于不存在非零、全非负的实系数 $\lambda_i$ 使 $\sum_i\lambda_iv_i=0$，可在上述二次域中精确判定。对下列各链，相邻顶点相同时定义为零；其余分支如下。

固定 $J=e_1$，边链为
\[
E(a,b)=\begin{cases}[a,b],&S(a,b),\\[a,aJ]+[aJ,b],&\text{否则}.
\end{cases}
\]
对 $Z=\sum_{\nu=1}^M c_\nu[v_{\nu0},\ldots,v_{\nu r}]$，定义
\[
t_* =\min\{t\in\{0,\ldots,3M\}:a(1,t,t^2,t^3)
\notin\operatorname{span}(v_{\nu0},\ldots,v_{\nu r})\ \forall\nu\},
\]
\[
C_a(Z)=\sum_\nu c_\nu[a(1,t_*,t_*^2,t_*^3),v_{\nu0},\ldots,v_{\nu r}],
\qquad C_a(0)=0.
\]
然后明确规定
\[
F(a,b,c)=\begin{cases}
[a,b,c],&S(a,b,c),\\
C_a(E(b,c)-E(a,c)+E(a,b)),&\text{否则},
\end{cases}
\]
\[
T(a,b,c,d)=\begin{cases}
[a,b,c,d],&S(a,b,c,d),\\
C_a(F(b,c,d)-F(a,c,d)+F(a,b,d)-F(a,b,c)),&\text{否则}.
\end{cases}
\]
每次候选搜索均有上界：一个待避真子空间在三次矩曲线上至多排除三个整数。故 $F$ 至多 6 项、19 个候选；$T$ 至多 24 项、73 个候选。被锥起的各单形是安全单形，锥顶点在其线性张成外，故新增单形仍然安全。降秩项保留在链中，只将其三维体积记为零。

把上述**由群元素唯一确定**的链展开为

\[
T(e_0,q_1,q_1q_2,q_1q_2q_3)
=\sum_{\nu=1}^{M\le24}c_\nu[v_{\nu0},v_{\nu1},v_{\nu2},v_{\nu3}].
\]

每项的边长是 $\ell_{\nu,ij}=\arccos\frac{v_{\nu i}\cdot v_{\nu j}}{\|v_{\nu i}\|\|v_{\nu j}\|}$，按式 (2) 的六棱次序排列；$D_\nu=\det(v_{\nu0},\ldots,v_{\nu3})$。所求的三个 E 型公式统一写为

\[
\boxed{\displaystyle
\omega_{E_m,k}(q_1,q_2,q_3)=
\exp\!\left[-\frac{ik}{\pi}
\sum_{\nu:D_\nu\ne0}c_\nu\operatorname{sgn}(D_\nu)
\mathcal B(\ell_{\nu1},\ldots,\ell_{\nu6})\right],
\quad q_j\in\mathcal Q_m,\quad m=6,7,8.}
\tag{10}
\]

式 (1)、(3)–(7)、(9) 和本节的有限分支定义了式 (10) 的每个符号。这是三类 E 群的一个完整、可不连续的几何代表；没有把未知函数“Vol”留给程序处理。它采用通用有限分支，**尚未化成类似 D 型八行表那样的专属短表**；如果“公式”特指那种短表，这一进一步化简仍未交付。

其闭性来自 $\partial T=\delta F$ 以及 $\int_{S^3}dV=2\pi^2$；左乘保持方向、秩和候选次序，所以五边形五项形成整数周期的闭三链。这里的证明与下文有限样本验证分开。

### 3.3 三组可以手算的 E 型相位

共同取 $g_1=e_1=i\sigma_1$、$g_2=e_3=i\sigma_3$，则 $p_2=e_2$。每行给出完整四元数 $g_3$：

| 群 | $g_3$ | $p_3=g_1g_2g_3$ | $\omega_{A,k}$ |
|---|---|---|---|
| $2T$ | $(1,1,-1,-1)/2$ | $(1,1,1,1)/2$ | $e^{-ik\pi/32}$ |
| $2O$ | $(1,1,0,0)/\sqrt2$ | $(0,0,1,1)/\sqrt2$ | $e^{-ik\pi/16}$ |
| $2I$ | $(\varphi/2,\varphi^{-1}/2,0,-1/2)$ | $(0,1/2,\varphi/2,\varphi^{-1}/2)$ | $\exp[-\frac{ik}{2}\arctan\frac1{4\varphi+1}]$ |

这三行不是从输出小数猜出来的：第一行将正交球面四面体（体积 $\pi^2/8$）从中心分为四个全等部分；第二行沿 $e_2e_3$ 中点平分该四面体。第三行是 $e_0$ 与 $S^2$ 三角形 $(e_1,e_2,p_3)$ 的球面连接；其面积为
\[
\Omega=2\arctan\frac{\varphi^{-1}}{3+\varphi}
=2\arctan\frac1{4\varphi+1}.
\]
积分 $\int_0^{\pi/2}\sin^2u\,du=\pi/4$ 给出 $V=\pi\Omega/4$，从而得到表中相位。第二、三行确实使用了超出共享 $2T$ 子群的元素。

## 4. 源码、数值证据和审计边界

本轮实际结果：26 项单元测试通过；24 个独立三重积分样本最大体积差为 $1.84742\times10^{-13}$；24 个独立边长/二面角比较最大差为 $2.118\times10^{-83}$。三个 E 群各三组五边形的最大差分别为 $1.143\times10^{-65}$、$9.108\times10^{-66}$、$8.464\times10^{-66}$。另外重跑 Q8 全部 4096 个四元组（最大差 $1.138\times10^{-86}$）和 E 型元素精确闭包。数值精度和输入详见报告，不将这些经验差值当作严格误差界。

- `src/su2_omega.py`：`edge_lengths_quaternions` 实现式 (2)；`oriented_volume_from_edges` 实现 (3)–(7)；`omega_quaternions_edge` 使用 A 负号及全局链。原正号、二面角接口保留。
- `scripts/ade_phase.py`：`phase_from_exact_quaternions` 默认走边长路线；`phase_details` 返回精确输入、累计顶点、分支、每项系数/取向/六边长/体积、总相位。`method='angle'` 用于独立交叉检查。
- `scripts/check_ade_phase.py`：每个群三组明确四元组；JSON 保存五边形的 **全部五个相位及各自计算项**，以及左右乘积与残差。另检验三组手算例、中心 $(-1,-1,-1)$ 相位为 $-1$、含单位元的归一化例。
- `SU2_ADE_phase_validation.json`：本轮实际执行的详细数值报告。它不是对 $24^4,48^4,120^4$ 个四元组的穷举，也不提供区间算术误差证书。

实际复现命令（先安装 `.[validation]`）：
```bash
python -m unittest discover -s tests -v
python -m scripts.su2_edge_crosscheck
python -m scripts.check_ade_phase
python -m scripts.check_ade_review
make adepdf
make completepdf
```

审计发现并修正：

1. `dic_quaternion` 的 $a^rb$ 坐标末分量错误写成 $+\sin(\pi r/n)$，与负叉积乘法不符；应为负号。它此前未被八分支表调用，故不能用表格通过来覆盖这个错误。
2. E 型入口及 D 型相位包装中的 `int(k)` 会把 $k=0.5$ 变为 0，已改为拒绝非整数；E 型同时按 $|k|$ 增加内部精度。
3. A 约定包装在默认 mpmath 精度下取倒数，可能把请求的 55 位结果降为约 16 位；已在内部高精度直接计算负号相位。
4. 验证依赖遗漏 SymPy；现已在安装项和 requirements 中声明。
5. ADE 文稿的 $\partial T$ 最后一项误写成 $F(a,b)$，已修正为 $F(a,b,c)$。旧文稿前后两处 E 型残差不一致，现统一指向逐项报告。
6. `normalized_numeric` 曾将符号 `s2/s5` 原样传给浮点转换，导致异常；现在先代入正平方根。

仍需区分的证明范围：E 型类阶 24、48、120 依赖几何代表与基本 CS/String 类的识别及 [Epa–Ganter 定理 1.1](https://arxiv.org/pdf/1605.09192)，不由几个相位小数推出；几何代表的逐点值也不必是对应阶的单位根。旧交接关于任意 $n$ 的 D 型八分支表，只给出“使用进位恒等式可约去”的概述，未展示全部符号化消去或对应几何链的逐分支推导；本轮不将其当作已经补齐的独立证明。源 A 原附件当前不在 checkout，本轮沿用仓库明确的 A 负号约定，不声称重新核验缺失附件。
