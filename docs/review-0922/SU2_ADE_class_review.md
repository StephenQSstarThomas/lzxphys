# SU(2) ADE 类识别独立审查（2026-09-22，中间版本）

> 本文件保留中间审查意见。Dic 任意 (n) 的表格闭性、边链比较和 E 型类识别
> 的最终限定已整合到 `SU2_ADE_HANDOFF_FINAL.md`；不要把本文件的“最小修正清单”
> 当作当前未修复项。

审查对象：docs/review-0922/SU2_ADE_progress.md、scripts/dic_formula.py、
scripts/check_dic_table.py 以及 SU2_anomaly_ADE_agent_prompt_zh(1).md。
本文不修改主文件，只判断 A 型符号、\(C_n\) 非恰当性、\(\mathrm{Dic}_n\) 八分支表
和 E 型类识别所需的逻辑条件。

## 总体结论

1. A 的 Eq.~(7.1) 负号是正确的，前提是相位定义为源笔记更新页中指定的
\(Z_{\rm before}/Z_{\rm after}\)，且 \(\operatorname{Vol}_A(C)=\int_CdV_{S^3}\)
使用同一个有向链。它不是把 level \(k\) 偷换成 \(-k\)。
2. \(C_n\) 的非恰当性论证是有效的，但必须明确：进位函数代表确实是当前
SU(2) 几何代表的限制，或用完整的边、面链比较证明二者同类。
3. scripts/dic_formula.py 的八分支有明确且一致的坐标约定，现有有理数检查通过；
不过 scripts/check_dic_table.py 只检查 \(n=2,3,4,5\)，这不是任意 \(n\) 的证明。
“来自半球构造”必须补成一条统一的 \(\mathrm{Dic}_n\) 边面三链构造，或把
任意 \(n\) 的闭性降级为尚未证明。另一个实质遗漏是：A 的 Eq.~(4.10)--(4.14)
只处理固定边换面，未覆盖源表与 E/F/T 之间可能不同的边链。
4. E 型周期结论可由 Epa--Ganter 的 Theorem~1.1、Proposition~4.1 可靠得到，
但必须显式写出“我们的 \(\omega_A\) 是 \(\mathrm{String}(3)\) 类的
\(\pm k\) 倍”这一识别假设及其归一化理由。E 型有限元素闭包与少量五边形数值
核查本身不能证明类阶。

## 1. A Eq.~(7.1) 的负号

源笔记 A 的公式链是自洽的：
\[
\omega_3^A=-\frac{1}{24\pi^2}\operatorname{Tr}(u^{-1}du)^3
          =-\frac{dV_{S^3}}{2\pi^2},
\qquad \int_{S^3}\omega_3^A=-1
\]
（A 的式 (1.1)、(2.12)），并且
\[
\mathcal A(A_0;\Lambda)=-\frac{k}{2\pi^2}\operatorname{Vol}_A(C),
\qquad
\omega_k^A=\exp\!\left(-\frac{ik}{\pi}\operatorname{Vol}_A(C)\right)
\]
（式 (2.14)、更新后的式 (7.1)）。因此
\[
\omega_k^A=(\omega_k^{\rm old})^{-1}
\]
只有在两者使用同一有向链、同一整数 \(k\)、而网络比值分别取源笔记
的 \(Z_{\rm before}/Z_{\rm after}\) 与旧仓库的相反 F-move 约定时才成立。
最小的共同约定表应至少写明这三项：链取向、初末网络比值、CS 三形式的符号。

A 的 \(U(1)\subset SU(2)\) 校准进一步确认了符号。对
\(h(\alpha)=e^{i\alpha\sigma_1}\) 和
\(N(\alpha,\beta)=\lfloor(\alpha+\beta)/2\pi\rfloor\)，A 式 (6.20)--(6.21) 给出
\[
\operatorname{Vol}_A(T(\alpha_1,\alpha_2,\alpha_3))
  =-\pi\alpha_1N(\alpha_2,\alpha_3),
\]
\[
\omega_k^A(h(\alpha_1),h(\alpha_2),h(\alpha_3))
 =e^{\,ik\alpha_1N(\alpha_2,\alpha_3)}.
\]
所以对 \(h=e^{2\pi i\sigma_1/n}\)，
\[
\omega_{k,C_n}(h^a,h^b,h^c)
 =\exp\!\left(\frac{2\pi ika}{n}
       \left\lfloor\frac{b+c}{n}\right\rfloor\right).
\]
这里正指数来自 \(\operatorname{Vol}_A<0\) 与 Eq.~(7.1) 的负号共同作用；
不能据此把 A 的 Eq.~(7.1) 改成正号。

相应地，若 \(I_4=-k\,\operatorname{Tr}(F^2)/(8\pi^2)\) 的局部 CS 代表
使用正体积 \(\eta_+=dV/(2\pi^2)\) 和正向旧 F-move，则其相位是旧仓库的
\(e^{+ikV/\pi}\)；A 的 \(Z_{\rm before}/Z_{\rm after}\) 是其逆代表。
最终报告必须保留“同一链时取逆”的限定，避免把网络方向转换误写成 level 转换。

## 2. \(C_n\) 的非恰当性

令
\[
z_n=\sum_{r=0}^{n-1}[h\mid h^r\mid h].
\]
在规范化 bar 复形中，边界两组项望远镜相消，故 \(\partial z_n=0\)。
对上面的进位函数代表，
\[
\prod_{r=0}^{n-1}\omega_{k,C_n}(h,h^r,h)
 =\exp(2\pi ik/n),
\]
因为只有 \(r=n-1\) 的进位
\(\lfloor(r+1)/n\rfloor\) 非零。若
\(\omega_1=\delta\beta\) 在整个离散化 SU(2) 上恰当，则限制满足
\[
\omega_1|_{C_n}=\delta(\beta|_{C_n}).
\]
任意二上链余边界在三循环上的乘积为 \(1\)，故 \(n=4\) 时
\[
P_4(\omega_1)=e^{2\pi i/4}=i\ne1
\]
立即矛盾。该论证不需要先分类 \(H^3(SU(2)^\delta,U(1))\)。

这里真正需要补出的同类性说明是：进位函数公式必须是当前几何 \(\omega_A\)
在 \(C_n\) 上的限制，或与它相差一个明确的二上链余边界。A 的式
(4.10)--(4.14) 在固定边链时可完成这一比较，但当前 E/F/T 规则的边
\(\gamma_{\rm EFT}\) 是径向短边（对反足点还使用固定 \(j=i\sigma_1\) 的分段边），
而 A 的 \(C_n\) 构造使用 \(\gamma_A(h(\alpha))(s)=h(s\alpha)\)。
当 \(\alpha>\pi\) 时二者一般不是同一边链；因此只引用换面三链还不够。

最小严格修正如下。取
\[
\partial A(g)=\gamma_A(g)-\gamma_{\rm EFT}(g),
\]
再选 \(B(g,h)\) 满足
\[
\partial B(g,h)=S_A(g,h)-S_{\rm EFT}(g,h)
-\bigl(gA(h)-A(gh)+A(g)\bigr).
\]
右端是二循环，由 \(H_2(S^3;\mathbb Z)=0\) 可填；随后
\[
\partial\bigl(C_A-C_{\rm EFT}-\delta_GB\bigr)=0.
\]
定义
\[
\beta(g,h)=\exp\!\left(2\pi ik\int_{B(g,h)}\eta_A\right)
\]
（若把 \(\eta_A=-dV/(2\pi^2)\) 吸收入 \(\eta_A\)，符号随之固定），即得
\[
\omega_A=\omega_{\rm EFT}\,\delta\beta
\]
或其等价的逆比值形式，取决于 A 链与 EFT 链的先后顺序。
这才把 \(P_4=i\) 证明为当前具体代表的限制不变量，而不是另一个已知
\(C_4\) cocycle 的不变量。

## 3. \(\mathrm{Dic}_n\) 八分支表

源表变量是
\[
g_i=V^{s_i}U^{r_i},\qquad
U=e^{i\pi\sigma_1/n},\quad V=i\sigma_2,\quad
r_i\in\{0,\ldots,2n-1\},\quad s_i\in\{0,1\}.
\]
这里 \(U\) 的阶为 \(2n\)，所以 \(\langle U\rangle=C_{2n}\)；这不要和
上一节的 \(C_n=\langle e^{2\pi i\sigma_1/n}\rangle\) 混淆。
表中的
\[
N_n(x,y)=\left\lfloor
\frac{[x]_{2n}+[y]_{2n}}{2n}\right\rfloor
\]
是 \(U\) 的 \(2n\) 进位。

脚本的 DicElement(n,r,epsilon) 表示 \(U^rV^\epsilon\)，不是表的
\(V^\epsilon U^r\)。当 \(\epsilon=1\) 时
\[
U^rV=VU^{-r},
\]
因此 table_coordinates 返回 \(([-r]_{2n},1)\) 是正确的；
这是八分支表与脚本输入标签之间不可省略的转换。
在这一转换下，表中
\[
\omega_A=\exp\!\left[-i\pi k\,
  \frac{\operatorname{Vol}_A}{\pi^2}\right]
\]
与 dic_phase 的实现一致。

check_dic_table.py 的精确有理核验检查
\[
\delta v_n\in2\mathbb Z,
\qquad
\delta v_n=v_n(b,c,d)-v_n(ab,c,d)+v_n(a,bc,d)
-v_n(a,b,cd)+v_n(a,b,c),
\]
这正是对所有整数 \(k\) 保证
\(\delta\omega_A=\exp(-i\pi k\delta v_n)=1\) 的充分条件。
它对 \(n=2,3,4,5\) 的全部四元组通过；我另以同一脚本逻辑扩展到
\(n=6,\ldots,10\)，也得到偶整数，但这仍是有限计算证据，不是任意 \(n\) 的证明。

当前文字称八分支表“按源笔记 A 的半球构造”，但 A 正文只显式推导了
最大环面 \(U(1)\) 的边、半球面和三链（式 (6.3)--(6.21)），没有在给定资料中
逐分支推导非交换 \(V^\epsilon U^r\) 的八个表达式。因此：

- 表式与所有已运行检查一致，不能指出一处数值公式错误；
- 有限穷举不能替代“对任意 \(n\)”的证明；
- 最小补充应二选一：给出八个 \(s_1s_2s_3\) 分支从统一边、面、三链
  \(\partial C_A=\delta_GS_A\) 的参数化推导；或逐分支用
  \[
  V^sU^r\,V^tU^q
  =V^{s+t\bmod2}U^{\,(-1)^t r+q+n st}
  \]
  （指数模 \(2n\)）及进位恒等式展开，证明每个 \(\delta v_n\) 是偶整数。

若采用第一种补充，闭性不再依赖穷举：统一链边界立即给出闭性，八个表式只是
该链积分的八种输入类型。若采用第二种补充，需把所有八种三输入与四输入位型
的余数/进位分支写入附录或提供可审计的符号化证明。在补充完成前，进度稿中的
“任意 \(\mathrm{Dic}_n\) 已精确闭性核验”应改为“任意输入由八分支表求值；
\(n=2,\ldots,5\) 已穷举核验，任意 \(n\) 的闭性需统一链或符号化证明”。

八分支短表与 E/F/T 使用不同面甚至不同边选择时，也必须使用上一节的
\(A(g),B(g,h)\) 比较。若实际证明两套边链相同，才可退回 A 的
式 (4.10)--(4.14)；当前材料没有给出该等同性，故不应只写“换面链产生
\(\delta\beta\)”。

## 4. E 型类识别与周期

令 \(\Gamma\subset S^3=\SU(2)=\Spin(3)\) 为有限子群，按离散群看待，
且作用在 \(S^3\) 上是自由的左乘作用。二元四面体 \(2T\)、二元八面体
\(2O\)、二元二十面体 \(2I\) 都满足这些条件，阶分别为 \(24,48,120\)。

对有限 \(\Gamma\)，实系数群上同调在正次数消失（平均收缩），所以系数列
\[
0\to\mathbb Z\to\mathbb R\to U(1)\to0
\]
给出 Bockstein 同构
\[
H^3_{\rm gp}(\Gamma,U(1))
\cong H^4(B\Gamma;\mathbb Z).
\]
有限自由球面群的周期上同调为
\[
H^4(B\Gamma;\mathbb Z)\cong\mathbb Z/|\Gamma|\mathbb Z.
\]
Epa--Ganter, arXiv:1605.09192v1, p.~12 的周期性陈述明确给出
正次数 \(4\) 的 \(\mu_{|\Gamma|}\)；这里 \(\mu_{|\Gamma|}\) 表示循环群
\(\mathbb Z/|\Gamma|\mathbb Z\)。

该文 Theorem~1.1（p.~2）说
\[
\operatorname{ord}\bigl(\mathrm{String}(3)|_\Gamma\bigr)=|\Gamma|
\]
对每个有限 \(\Gamma\subset S^3\) 成立。Proposition~4.1（pp.~12--13）进一步给出
\[
\mathrm{String}(3)|_\Gamma
\simeq\Gamma_{\rm uni}
[\mu_{|\Gamma|}\hookrightarrow U(1)].
\]
其证明使用纤维化
\[
S^3/\Gamma\longrightarrow B\Gamma\longrightarrow
BS^3\simeq\mathbb H P^\infty,
\]
Leray--Serre 谱序列中唯一非平凡 \(d_4\) 为乘以 \(|\Gamma|\)。
Remark~4.3（p.~13）特别核对了循环子群和最大环面 \(S^1\subset S^3\) 的
限制及其符号。

把上述定理用于本文需要显式添加以下类识别引理：
\[
[\omega_A|_\Gamma]=\varepsilon k
[\mathrm{String}(3)|_\Gamma],\qquad \varepsilon=\pm1.
\]
这不是元素闭包或少量五边形数值自动给出的。它依赖普通带联络 CS
特征类的自然性，以及 A 的归一化
\[
\omega_3^A=-\frac{dV}{2\pi^2},\qquad
\int_{S^3}\omega_3^A=-1,
\]
说明 A 类与正向体积生成元只差逆元；A 的 Eq.~(2.14)/(7.1) 再固定网络比值的
\(\varepsilon\)。若前述参考三骨架截面与边、面链比较已证明 \(\omega_A\) 是
\(-k\widehat c_2\) 的平坦限制，则此引理成立。否则 Epa--Ganter 只证明
String 类的阶，不能单独证明当前 E 型代码输出的类阶。

一旦该识别成立，令 \(x_\Gamma\) 为 String(3) 限制的生成元，则
\[
[\omega_k|_\Gamma]=\varepsilon kx_\Gamma,\qquad
\operatorname{ord}[\omega_k|_\Gamma]
=\frac{|\Gamma|}{\gcd(k,|\Gamma|)},
\]
\[
[\omega_k|_\Gamma]=0\quad\Longleftrightarrow\quad
|\Gamma|\mid k.
\]
取逆只改变生成元符号，不改变阶与 \(k\) 周期。这里的类结论不等于未重定相
代表逐点满足 \(\omega_k(g,h,l)^{|\Gamma|}=1\)；它只表示
上同调类的 \(|\Gamma|\) 倍为零。

## 5. 附带实现问题

表格代码本身的坐标转换、相位负号和 \(\delta v\in2\mathbb Z\) 检查逻辑一致。
另有一个不在本轮表格调用路径中的文档/实现不一致：
scripts/ade_subgroups.py 的 dic_quaternion 文档写成 \(a^r b^\epsilon\)，
但 \(\epsilon=1\) 返回
\((0,0,\cos(\pi r/n),+\sin(\pi r/n))\)，这是 \(ba^r=VU^r\) 的坐标；
按当前 \(+i\sigma\) 的负叉积，直接计算 \(a^rb\) 应为
\((0,0,\cos(\pi r/n),-\sin(\pi r/n))\)。该函数目前未被
dic_formula.py 调用，故不改变八分支核验；若作为公开嵌入入口，应修正文档
或返回值，并保留 table_coordinates 的 \(r\mapsto-r\) 转换。

## 6. 最小修正清单

1. A 符号段加入同一有向链及同一 \(Z_{\rm before}/Z_{\rm after}\) 比值的限定，
   并引用 A 式 (1.1)、(2.12)、(2.14)、(6.20)--(6.21)、(7.1)。
2. \(C_n\) 段加入边链比较；否则只能证明 A 半球代表非恰当，不能自动转给 E/F/T。
3. Dic 表保留现有有限穷举证据，但把任意 \(n\) 的闭性改由统一半球链推导或
   可审计的进位恒等式证明承担。
4. “与 E/F/T 同类”加入 \(A(g)\) 和 \(B(g,h)\) 的边、面比较，除非先证明两套
   \(\gamma\) 相同。
5. E 型段先写 \([\omega_A|_\Gamma]=\pm k[\mathrm{String}(3)|_\Gamma]\)，
   再引用 Epa--Ganter Theorem~1.1、Proposition~4.1、Remark~4.3 推出周期。

参考：N. Epa and N. Ganter, “Platonic and alternating 2-groups,”
arXiv:1605.09192v1 (2016), Theorem~1.1, Proposition~4.1, Remark~4.3,
\url{https://arxiv.org/abs/1605.09192}.
