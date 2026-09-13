# SU(2) 体积三余循环的 CS 类识别与符号审查

日期：2026-09-13。限域任务：证明指定几何 bar 链三余循环实现指定反常多项式的二级 Chern–Simons 类，修复全局纯规范及 C8/C11 推理，给出换面余边界与独立整数 level 校准。本文不计算整个离散化群上同调，不处理球面四面体的 Li₂ 求值。

## 1. 可直接采用的结论

令 \(G=SU(2)\)，使用反厄米联络、二维基本表示普通迹，以及题设外法向球面取向：

\[
P_k(F)=I_4=-\frac{k}{8\pi^2}\operatorname{Tr}(F\wedge F),\qquad
I_3(A)=-\frac{k}{8\pi^2}\operatorname{Tr}(A\wedge dA+\tfrac23A^3),
\]
\[
\eta=\frac{\operatorname{Tr}(h^{-1}dh)^3}{24\pi^2}
=\frac{dV_{S^3}}{2\pi^2},\qquad \int_{S^3}\eta=1.
\]

固定左平移等变的边、面、三维积分链，满足

\[
\partial\gamma(g)=[g]-[e],\qquad
\partial S(g,h)=g\gamma(h)-\gamma(gh)+\gamma(g),
\tag{1}
\]
\[
\partial C(g,h,l)=gS(h,l)-S(gh,l)+S(g,hl)-S(g,h).
\tag{2}
\]

则

\[
\boxed{\omega_k(g,h,l)=\exp\left(2\pi i k\int_{C(g,h,l)}\eta\right)}
\tag{3}
\]

是普通带联络 SU(2) Chern–Simons 微分特征类 \(-k\widehat c_2\) 的平坦限制的一个具体、可不连续的 bar 代表。这里的“具体”包括边和面的选定；改变这些选择一般改变代表，但不改变类。第 4 节给出逐胞腔的全局识别证明，不以局部三次展开推断该结论。

在本文标准 Bockstein 约定下，另外有

\[
\boxed{c_2=+\frac{\operatorname{Tr}(F^2)}{8\pi^2},\quad
u_k=-kc_2,\quad
\mathcal B([\omega_k])=+k\,\iota^*c_2,}
\tag{4}
\]

其中将 \(U(1)\) 通过 \(t\mapsto e^{2\pi it}\) 识别为 \(\mathbb R/\mathbb Z\)，
\(\mathcal B([t\bmod\mathbb Z])=[\delta t]\)，而
\(\iota:BG^\delta\to BG\)。式 (4) 中最后一个正号不是 \(u_k\) 的符号；原因是标准平坦微分特征满足 \(\mathcal B(a)=-c(a)\)，见第 3 节。

设 \(h=\exp(2\pi i\sigma_1/n)\)、\(n\ge2\)。对任意实现 (3) 的规范化全局代表，独立的有限子群校准为

\[
\boxed{\prod_{r=0}^{n-1}\omega_k(h,h^r,h)
=\exp(-2\pi i k/n).}
\tag{5}
\]

这既检验整数 level，也检验全局分支；若仅将所有降维输入设为 1，就会在 \(n\nmid k\) 时违反 (5)。

## 2. 有限规范变换与 F-move 的方向

置 \(A^u=u^{-1}Au+u^{-1}du\)。直接展开 Chern–Simons 三形式，并用 Maurer–Cartan 方程，得到

\[
\operatorname{Tr}(A^u dA^u+\tfrac23(A^u)^3)
=\operatorname{Tr}(A dA+\tfrac23A^3)
-d\operatorname{Tr}(du\,u^{-1}\wedge A)
-\tfrac13\operatorname{Tr}(u^{-1}du)^3.
\]

因而在当前普通迹及 level 约定下，精确的有限变换公式是

\[
\boxed{I_3(A^u)-I_3(A)
=dQ_2(A,u)+k\,u^*\eta,\qquad
Q_2(A,u)=\frac{k}{8\pi^2}\operatorname{Tr}(du\,u^{-1}\wedge A).}
\tag{6}
\]

这同时固定 Wess–Zumino 项和边界项的符号。特别地，在一个可缩且联络平坦的三维单纯形内，可选平行标架使初始联络为 0；以局部 developing map \(f\) 换到指定边界标架后，

\[
A=f^{-1}df,\qquad I_3(A)=k f^*\eta.
\tag{7}
\]

二维反常的有限规范变换通过三维插值上的指数化 CS 求值；其相对边界标架决定该值。对于只改变局部融合次序的 F-move，取

\[
B_L=S(g,h)+S(gh,l),\qquad
B_R=gS(h,l)+S(g,hl).
\]

它们有同一个一维边界，且 (2) 正是

\[
\partial C=B_R-B_L.
\tag{8}
\]

规定正向 F-move 为
\((D_gD_h)D_l\to D_g(D_hD_l)\)，并规定其标量是从 \(B_L\) 到 \(B_R\) 的正向 CS 插值相位，则 (7)–(8) 给出 (3)。若采用反向 F-move，或把“边界反常”定义为抵消当前正向 bulk 变化的因子，所写标量相应取逆；不能只改变这项约定而保持 (3) 的符号不变。

为什么这里没有遗漏 (6) 的 \(Q_2\)？单个局部插值在平行标架下从 \(A=0\) 开始，故该项为零。把相邻单纯形的结果拼起来时，边界标架的选择已由共享的边面数据固定；这正是 CS 的相对边界平凡化。更换它们产生第 6 节的结点重定相 \(\delta\beta\)。在一般背景标架直接计算，也必须保留 (6) 的边界项，结果与这一选择一致。

原论文 v2 的式 (7)–(10) 及附录 E，尤其 \(E6\)–\(E10\)，采用此有限规范变换／三维插值路线。与其正文规范化相配，应取其双线性型 \(\mathrm{tr}(XY)=-\operatorname{Tr}(XY)/(4\pi^2)\)。本文的全局步骤用下节微分特征类补足。[Jia 等，原论文 v2](https://arxiv.org/pdf/2510.14722v2)

## 3. 带联络的二级类，以及 C8/C11 的正确替代

按总陈类定义展开 \(\det(1+iF/(2\pi))\)，对 \(\operatorname{Tr}F=0\) 有

\[
c_2(F)=-\tfrac12\operatorname{Tr}\left(\frac{iF}{2\pi}\right)^2
=\frac{\operatorname{Tr}F^2}{8\pi^2}.
\tag{9}
\]

因此 (3) 应识别到积分特征类 \(u_k=-kc_2\)，而不是未经换号的 \(+kc_2\)。

用现代次数记法，取一个 degree-4 微分特征

\[
\widehat u_k(P,A):Z_3(M;\mathbb Z)\longrightarrow\mathbb R/\mathbb Z,
\]

其定义性质为

\[
\widehat u_k(P,A)(\partial b)=\int_bP_k(F)\pmod{\mathbb Z},
\quad \operatorname{curv}(\widehat u_k)=P_k(F),
\quad c(\widehat u_k)=u_k(P).
\tag{10}
\]

这里还要求对保持联络的束映射自然。Cheeger–Simons 定理 2.2 给出这种特征的存在和自然唯一性；命题 2.8 给出在总空间上用 CS 三形式表示、再经截面拉回的公式。原文次数为“作用于 3-循环的 degree-3 character”，与此处现代 degree-4 只是记号相差一。[Cheeger–Simons，定理 2.2、命题 2.8，印刷页 59–61](https://math.mit.edu/juvitop/pastseminars/notes_2019_Fall/cheeger-simons.pdf)

若 \(F=0\)，(10) 对边界为零，因此得到

\[
a_k(P,A)\in H^3(M;\mathbb R/\mathbb Z).
\tag{11}
\]

所有平坦束的自然性产生通用类

\[
a_k^\delta\in H^3(BG^\delta;\mathbb R/\mathbb Z),\qquad
a_k(P,A)=\phi_\delta^*a_k^\delta.
\tag{12}
\]

正确的数据是 \((P,A)\)，或等价的平坦束分类映射 \(\phi_\delta\)，不能仅保留 \(u_k(P)\)。这是 C11 应采用的形式。

为检查 Bockstein 符号，取 \(\widehat u_k\) 的实上链提升 \(t\)。微分特征条件写成

\[
\delta t=P_k(F)-c,\qquad c\in Z^4(M;\mathbb Z),\quad[c]=u_k(P).
\]

平坦时 \(\delta t=-c\)，所以

\[
\boxed{\mathcal B(a_k(P,A))=-u_k(P).}
\tag{13}
\]

该负号可直接核对 Cheeger–Simons 推论 1.2(1) 与平坦束公式 (8.2)、(8.5)。采用 (10) 的正向 holonomy 与通常 Bockstein 时，不能原样沿用论文 \(C8\) 的正号又把其 \([Tf]\) 直接认作当前 CS 相位；若使用相反 Bockstein 或相反二级类，必须明确宣布。[Cheeger–Simons，印刷页 54、70–71](https://math.mit.edu/juvitop/pastseminars/notes_2019_Fall/cheeger-simons.pdf)

仅在普通系数长正合列中给出一个 Bockstein 提升不能选定 (12)：两个提升仍可相差
\(\operatorname{im}(H^3(BG^\delta;\mathbb R)\to H^3(BG^\delta;\mathbb R/\mathbb Z))\)。本文选定的是由 (10) 和联络自然性确定的提升。

## 4. 从通用平坦束到指定 bar 链的逐胞腔证明

### 4.1 通用束和三骨架截面

采用齐次 bar 模型：\(EG^\delta\) 的 \(p\)-单纯形写成
\([a_0,\ldots,a_p]\)，左作用为同时左乘各顶点。商 \(BG^\delta\) 的非齐次坐标为
\(g_i=a_{i-1}^{-1}a_i\)。于是规范化三单纯形的累计顶点为

\[
[e,g,gh,ghl].
\tag{14}
\]

考虑主右 \(G\)-束

\[
\mathcal P=(EG^\delta\times G)/G^\delta,
\qquad a\cdot(x,v)=(ax,av),\qquad [x,v]\cdot r=[x,vr].
\tag{15}
\]

这是沿 \(G^\delta\to G\) 扩张结构群的通用平坦束。为明确左右作用，\(EG^\delta\) 本身的主右作用为 \(x\cdot a=a^{-1}x\)；映射 \(x\mapsto[x,e]\) 是右等变的，因而 (15) 的普通分类映射确为 \(\iota\)。

这里的联络和截面按各个光滑单纯形及其相容面来读。\(BG^\delta\) 本身不是有限维光滑流形；通用微分特征的说法是 (12) 所定义的自然平坦类。以下积分可等价地在有限三角剖分的循环上进行，或拉回到带平坦束的光滑流形上计算，再由自然性作通用识别。

乘积上联络取纤维 Maurer–Cartan 形式

\[
\Theta=v^{-1}dv.
\]

同时左乘常数 \(a\) 不改变 \(\Theta\)，故它下降到 \(\mathcal P\)，且 \(d\Theta+\Theta^2=0\)。

先取能够由逐单纯形映射实现的一套边面数据。构造左等变分段光滑映射

\[
f:(EG^\delta)^{(3)}\longrightarrow G,\qquad f([a])=a.
\tag{16}
\]

构造可逐轨道完成：连通性给出边；\(\pi_1(SU(2))=0\) 使每个边界圆延拓到二维面；\(\pi_2(SU(2))=0\) 使每个面边界球延拓到三单纯形。归一化退化单纯形取退化映射。等变性通过左平移定义，无需对群标签具有全局连续依赖，因为此处是 \(G^\delta\) 的 bar 模型。

把各单纯形基本链经 \(f\) 推前，得到 \(\gamma,S,C\)，其边界逐字满足 (1)–(2)。映射 (16) 给出三骨架上的截面

\[
s([x])=[x,f(x)].
\tag{17}
\]

这个定义与提升 \(x\) 的选择无关：\(f(ax)=af(x)\)。在该截面中

\[
A_s=s^*\Theta=f^{-1}df.
\tag{18}
\]

式 (18) 的右端在每个提升上书写，但由于左平移不变性，它拼成基底三骨架的联络形式。这里没有声称存在基底上的单值 \(f\)。

### 4.2 每个三胞腔上的值就是体积值

CS 在总束空间上的拉回公式与 (18) 给出：对三骨架中的任意整系数三循环
\(z=\sum_\sigma n_\sigma\sigma\)，

\[
\begin{aligned}
\widehat u_k(\mathcal P,\Theta)(z)
&=\int_z I_3(A_s)\pmod{\mathbb Z}\\
&=k\sum_\sigma n_\sigma\int_{f_*\sigma}\eta\pmod{\mathbb Z}\\
&=\sum_\sigma n_\sigma\,
\frac{\log\omega_k(g_{\sigma,1},g_{\sigma,2},g_{\sigma,3})}{2\pi i}
\pmod{\mathbb Z}.
\end{aligned}
\tag{19}
\]

最后一行中的“log”仅表示 \(\mathbb R/\mathbb Z\) 值，不需要选复对数分支。等式对所有三循环成立；它直接识别通用平坦 CS holonomy 与几何上链，不只识别它们的导数、局部展开或 Bockstein。

第 5 节说明右端对四维边界为零，故下降到 \(H_3(BG^\delta;\mathbb Z)\)。每个三维同调类都能表示在三骨架中。又因为 \(\mathbb R/\mathbb Z\) 是可除群，通用系数定理的 Ext 项为零，三上同调类由其对三循环的取值完全决定。因此

\[
\boxed{[\omega_k]=a_k^\delta.}
\tag{20}
\]

对任意光滑流形上的平坦束，把其分类映射拉回上述构造；在三角剖分的每个单纯形内取平行标架，累计群元给出 (14)，得到同一个逐四面体公式。有限规范变换的局部 F-move 相位因而就是通用类 (20) 在所选结点平凡化中的代表。

### 4.3 任意既定积分链仍得到同一个类

若题设的 \(S,C\) 是有限分段光滑积分链而非单个盘／球的推前，不必假装它们必能原样实现为 (16)。先以上述映射构造作为参考，再应用第 6 节的链比较公式；\(H_1(S^3)=H_2(S^3)=0\) 保证比较链存在。它显示任意满足 (1)–(2) 的预定数据只相差显式 \(\delta\beta\)，所以仍代表 (20)。这也处理以多个球面单纯形相加实现面、三链的全局算法。

几何链的顶点始终是 \(SU(2)\) 中的实际单位四元数。若实现暂存精确齐次有理向量，须在“定义链”时将每一顶点取正范数单位化；未经单位化的正缩放向量只是同一球面点的坐标数据，不是不同奇异链的顶点。本文所有边界等式都在实际球面链中成立。

## 5. 五边形、归一化及全局纯规范漏洞

在取值于群上积分链的 bar 复形中，记左作用余边界为 \(\delta_G\)。式 (1)–(2) 是
\(\partial S=\delta_G\gamma\)、\(\partial C=\delta_GS\)。对四个群元有

\[
Z=g_1C(g_2,g_3,g_4)-C(g_1g_2,g_3,g_4)
+C(g_1,g_2g_3,g_4)-C(g_1,g_2,g_3g_4)+C(g_1,g_2,g_3).
\]

边界算子与左平移可交换，故

\[
\partial Z=\delta_G\partial C=\delta_G^2S=0.
\]

这是真正共享、带符号的面相消。由于 \(\eta\) 左不变且具有整数周期，

\[
(\delta\omega_k)(g_1,g_2,g_3,g_4)
=\exp\left(2\pi i k\int_Z\eta\right)=1.
\tag{21}
\]

这里 \(Z\) 不一定是边界，也无需虚构把它填进四维球的映射；其整数周期已经足够。\(\pi_3(SU(2))\) 恰恰意味着三骨架截面一般不能延到四骨架，其四胞腔阻碍与这一整数相联系。

取 \(\gamma(e)=0\)、\(S(e,g)=S(g,e)=0\)。只要一个输入为单位元，(2) 的右端为零，此时 \(C\) 是三循环，指数就是 1；也可以直接选 \(C=0\)。因此 \(\omega\) 规范化。

固定面时，两个填充的差也是三循环，所以更换三维填充不改变指数。这一陈述不能推广成“更换面也不改变代表”。

关于原论文式 (10) 前的全局论证：SU(2) 束在三维基底上平凡，只保证全局截面和全局联络一形式，不保证存在基底上的 \(h\) 使 \(A=h^{-1}dh\)。例如在 \(S^1\times S^2\) 的平凡束上，\(A=i\lambda\sigma_3\,d\theta\) 平坦，但一般有非平凡 holonomy。本文只在可缩胞腔上使用纯规范；整体上使用 (15)–(18) 的等变 developing map 与 (10)，因此不需要错误的全局推论。

## 6. 更换面或边：显式 \(\beta\) 与符号

### 6.1 边固定，只换面

令 \(S'\) 和 \(S\) 使用同一套 \(\gamma\)。则
\(\partial(S'-S)=0\)，由 \(H_2(S^3;\mathbb Z)=0\)，可逐对群元选三链 \(B(g,h)\) 满足

\[
\partial B(g,h)=S'(g,h)-S(g,h).
\tag{22}
\]

设 \(C'\) 对应 \(S'\)。直接计算

\[
\partial\big(C'-C-\delta_GB\big)=0.
\]

因此

\[
\boxed{\beta(g,h)=\exp\left(2\pi i k\int_{B(g,h)}\eta\right),\qquad
\omega'_k=\omega_k\,\delta\beta,}
\tag{23}
\]
\[
(\delta\beta)(g,h,l)
=\frac{\beta(h,l)\beta(g,hl)}{\beta(gh,l)\beta(g,h)}.
\]

这里 \(\beta\) 前面的符号是正号，因为 (22) 取的是 \(S'-S\)。更换 \(B\) 的填充只添三循环，所以 \(\beta\) 本身也无填充歧义。

### 6.2 边也改变

取二链 \(A(g)\) 满足

\[
\partial A(g)=\gamma'(g)-\gamma(g)
\]

（由 \(H_1(S^3)=0\) 保证存在），再取三链

\[
\partial B(g,h)=S'(g,h)-S(g,h)
-\big(gA(h)-A(gh)+A(g)\big).
\tag{24}
\]

右端为二循环，故仍有解。应用 \(\delta_G^2=0\)，同样得
\(\partial(C'-C-\delta_GB)=0\)，从而 (23) 不变。必要时对单位元输入选零比较链以保留规范化。

在缺陷语言中，(23) 是三价结点选择的重定相。它不会改变任何闭合三循环的 CS holonomy，尤其不会改变 (5)。

## 7. 为何没有悄悄加上另一个反常类

这一点有两层，不能只由 \(\delta\omega=1\) 得出。

1. **普通带联络 SU(2) 反常的输入已唯一。** \(BSU(2)\simeq\mathbb H P^\infty\) 在相关次数只有 0、4 维胞腔，因此 \(H^4(BSU(2);\mathbb Z)=\mathbb Z c_2\)、\(H^3(BSU(2);\mathbb R/\mathbb Z)=0\)，且四次整系数到实系数的映射单射。指定 \(P_k\) 就指定了积分提升 \(-kc_2\)，没有遗漏同一实多项式下的扭结整类；指定联络自然性后，也没有另一个通用平坦微分特征可加。等价地可直接使用第 3 节所引自然唯一性定理。

2. **具体 bar 代表已在全部三循环上匹配。** 式 (19) 是通用 CS 类与体积上链对任意三循环的逐项等式；式 (23)–(24) 又覆盖全部允许的几何选择。因此，即使 \(H^3(BG^\delta;U(1))\) 还含有其他类，它们也没有进入本构造。无需分类该群，更不能只凭相同小角三次项或相同 Bockstein 排除它们。

如果只给平坦背景上的局部密度，并允许任意离散群上同调附加项，则该信息确实不足以确定全局类。本任务使用的是普通 SU(2) 带联络反常内流，再限制到平坦背景；这项物理输入正是上述唯一性适用的条件。

Baez–Lauda §8.5、定理 55 及印刷页 68–69 用带联络二级类建立整数族进入群三上同调的映射；其作用不是给 \(C8\) 的任意 Bockstein 提升补一句“规范”。本文的 (19) 补充到指定链代表层次。[Baez–Lauda v3](https://arxiv.org/pdf/math/0307200v3)

## 8. 有限环面子群：两个独立的符号校准

### 8.1 标准 bar 循环与 Bockstein 计算

将 \(C_n\) 元素写成 \(h^a\)，\(0\le a<n\)，记 \([a+b]_n\) 为标准余数，

\[
q(a,b)=\frac{a+b-[a+b]_n}{n}=\left\lfloor\frac{a+b}{n}\right\rfloor.
\]

整数二余循环 \(q\) 表示 \(x\in H^2(BC_n;\mathbb Z)\)。基本 SU(2) 表示限制为两个互逆权的线表示，故

\[
c_2|_{BC_n}=-x^2,\qquad u_k|_{BC_n}=+kx^2.
\tag{25}
\]

可能把权一线的 \(c_1\) 定义成 \(-x\) 的习惯不影响平方及 (25)。定义实三上链

\[
t_m(a,b,c)=\frac{m a}{n}q(b,c).
\]

用 \(\delta q=0\) 及 \(a+b-[a+b]_n=nq(a,b)\) 逐项相消可得

\[
\delta t_m(a,b,c,d)=m q(a,b)q(c,d).
\tag{26}
\]

故 \(\exp(2\pi it_m)\) 的标准 Bockstein 为 \(m x^2\)。有限群的正次数实系数上同调可通过平均收缩，因而此处 Bockstein
\(H^3(BC_n;\mathbb R/\mathbb Z)\to H^4(BC_n;\mathbb Z)\)
是同构。由 (13)、(25)，当前 CS 类对应 \(m=-k\pmod n\)。

规范化 bar 三链

\[
z_n=\sum_{r=0}^{n-1}[h\mid h^r\mid h]
\tag{27}
\]

是循环：应用标准边界后，\([h^r\mid h]-[h^{r+1}\mid h]\) 和
\([h\mid h^{r+1}]-[h\mid h^r]\) 分别望远镜相消。规范化复形中的单位元项为零。又

\[
\langle t_m,z_n\rangle
=\frac{m}{n}\sum_{r=0}^{n-1}q(r,1)=\frac{m}{n}.
\]

这证明 (5)，并表明余边界对该乘积没有影响。

### 8.2 直接扫半球，不借 Bockstein 决定负号

下面在 \(C_n\) 上选择一套允许的链，直接算出同一个代表；第 6 节保证所得循环不变量与其他一致选择相同。

令 \(E\) 是最大环面上从 \(e\) 正向到 \(h\) 的短弧，

\[
\gamma(h^a)=\sum_{j=0}^{a-1}h^jE,\qquad
T=\sum_{j=0}^{n-1}h^jE.
\]

这是长弧链约定，\(\gamma(e)=0\)；\(T\) 为整条正向大圆，且作为细分奇异链满足 \(h^aT=T\)。取

\[
H=\{(x_0,x_1,x_2,0)\in S^3:x_2\ge0\}
\]

为有向半球，使 \(\partial H=T\)，则可选

\[
S(h^a,h^b)=q(a,b)H.
\tag{28}
\]

令 \(K_a\) 为从 \(H\) 扫到 \(h^aH\) 的三链，并使

\[
\partial K_a=h^aH-H.
\tag{29}
\]

严格的奇异链构造如下：对同伦 \(L_{\exp(it\sigma_1)}\)，\(0\le t\le2\pi a/n\)，用棱柱算子 \(P_a\) 取
\(\partial P_aH=h^aH-H-P_aT\)。因 \(h^aT=T\)，\(P_aT\) 是取值于该大圆的二循环。由 \(H_2(S^1)=0\) 取同在大圆中的三链 \(R_a\)，使 \(\partial R_a=P_aT\)，再置 \(K_a=P_aH+R_a\)。\(R_a\) 对三形式的积分为零，因而它只保证严格链边界，不改变以下体积。

取半球参数 \(0\le\rho\le\pi/2\)、\(0\le\varphi<2\pi\)，半球取向为 \(d\varphi\wedge d\rho\)，其边界正是正向 \(T\)。左乘 \(e^{it\sigma_1}\) 后的参数化为

\[
x(t,\varphi,\rho)=
\big(\cos\rho\cos(\varphi+t),
\cos\rho\sin(\varphi+t),
\sin\rho\cos t,-\sin\rho\sin t\big).
\tag{30}
\]

此处法平面 \((x_2,x_3)\) 旋转角为 \(-t\)，正是 \(+i\sigma\) 乘法约定的负叉积。以 \(dt\wedge d\varphi\wedge d\rho\) 取向，直接算环境行列式得到

\[
x^*\eta=-\frac{\cos\rho\sin\rho}{2\pi^2}
dt\wedge d\varphi\wedge d\rho,
\qquad
\int_{K_a}\eta=-\frac{a}{n}.
\tag{31}
\]

由 \(\delta q=0\)，(28) 的 bar 边界恰为
\(q(b,c)(h^aH-H)\)，所以可选

\[
C(h^a,h^b,h^c)=q(b,c)K_a.
\]

将 (31) 代入 (3) 得到完全显式的有限限制代表

\[
\boxed{\omega_k^{C_n}(h^a,h^b,h^c)
=\exp\left(-\frac{2\pi i k a}{n}
\left\lfloor\frac{b+c}{n}\right\rfloor\right).}
\tag{32}
\]

于是 (5) 的负号可在几何中直接核验。正向商 \(C_n\backslash S^3\) 上 \(\int\eta=+1/n\) 并不与 (5) 冲突：将 \(\sigma_1\) 共轭到 \(\sigma_3\) 后，左作用在第一列的复坐标中具有权 \((1,-1)\)，通常记作 \(L(n,-1)\)。商束的 developing map 取球面上的恒等映射，故正向基本类的当前 CS 值是 \(+1/n\)；与 (32) 在 \(z_n\) 上的 \(-1/n\) 比较，得到分类映射送出的基本类是 \(z_n\) 的负类。直接把商球面的正向积分当成 (27) 的取值会丢掉该取向差。式 (30)–(31) 无需依赖 lens 空间命名习惯。

两个整数 level 若给出相同全局类，则对所有 \(n\) 应有相同的 (5)，从而每个 \(n\) 都整除其差，故二者相等。这在当前整数族内部给出独立单射性校准，不涉及 \(H^3(G^\delta;U(1))\) 的其余部分。

## 9. 交付状态

已给出：有限规范变换的精确边界项及符号；F-move 方向；通用平坦束／三骨架截面的逐胞腔类识别；任意允许积分链与参考截面的比较；五边形及归一化证明；换边换面时显式 \(\beta\)；普通 \(c_2\)、Bockstein 及整数 level 的双重独立符号校准。

辅助核查：用精确有理数穷举验证了 \(n=2,\ldots,9\) 的式 (26) 与循环配对式；用符号行列式验证了 (30) 的雅可比为 \(-\cos\rho\sin\rho\)。这些核查只检查代数及符号，证明仍由上述推导承担。

以上结论是链与微分特征层面的全局证明，不以数值五边形或小角展开为证据。任意输入的具体几何链如何有限构造、其体积如何化成指定特殊函数闭式，属于其他交付；本文不把这些求值工作视为已经完成。
