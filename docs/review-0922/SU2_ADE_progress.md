# 2026-09-22 SU(2) ADE 反常推进审计

本文件是对 `SU2_anomaly_ADE_agent_prompt_zh(1).md` 的第一轮执行记录。范围集中于三件事：新资料对已有 SU(2) 公式的约定修正、闭性与非恰当性的可验证证据、以及二元二面体群 (2D_n=mathrm{Dic}_n) 的精确标签和求值入口。E 型群的元素集合和闭包也已实现；全部 ADE 的上同调周期先用已核验的有限子群定理识别，尚未打印每个 E 型群的逐三元组相位表。

## 1. 新资料的关键修正：整体符号必须显式转换

新附件 A 的第 10 页更新式 (7.1) 明确采用

\[
 \omega_k^{A}(g_1,g_2,g_3)
 =\exp\!\left(-\frac{ik}{\pi}\operatorname{Vol}_{A}(C)
 \right),
 \qquad \operatorname{Vol}_{A}(C)=\int_C dV_{S^3}.
\]

同一附件式 (1.1)、(2.12) 的归一化是

\[
 \omega_3^A=-\frac{dV_{S^3}}{2\pi^2},
 \qquad \int_{S^3}\omega_3^A=-1.
\]

参考论文 E 的式 (10) 也给出左不变三形式的 CS 取值；它与 A 的符号可通过 `tr=Tr/(4π²)` 的约定核对。旧仓库文稿和 `src/su2_omega.py` 使用的是

\[
 \eta_+=+\frac{dV_{S^3}}{2\pi^2},qquad
 \omega_k^{\rm old}=\exp\!\left(+\frac{ik}{\pi}V_{\rm old}\right).
\]

因此本轮采用的明确换算由 `src/su2_omega.py` 的新入口
`omega_quaternions_source` 实现：

\[
 \boxed{\omega_k^A=(\omega_k^{\rm old})^{-1}}
\]

在两者使用同一有向几何链时成立。这里不是悄悄修改 level：它是“源笔记 A 的 (Z_{\rm before}/Z_{\rm after}) 约定”与旧仓库正向相位的代表转换。A 的 (U(1)) 校准式 (6.21) 给出

\[
 \omega_k^A(h(\alpha_1),h(\alpha_2),h(\alpha_3))
 =\exp\!\left[
 \frac{ik\alpha_1}{2\pi}
 \big(\alpha_2+\alpha_3-[\alpha_2+\alpha_3]_{2\pi}\big)
 \right].
\]

这与 (\operatorname{Vol}_{A}=-\pi\alpha_1N(\alpha_2,\alpha_3)) 和式 (7.1) 一致。旧仓库的 (C_n) 检查得到逆相位；本轮结果统一以 A 约定报告。

## 2. closed 性质：证明与实际 Q8 检查

沿用已经审查过的单位顶点规范化奇异链。对边、面、三维填充分别记 (E,F,T)，它们满足

\[
\begin{aligned}
 \partial E(a,b)&=[b]-[a],\\
 \partial F(a,b,c)&=E(b,c)-E(a,c)+E(a,b),\\
 \partial T(a,b,c,d)&=F(b,c,d)-F(a,c,d)+F(a,b,d)-F(a,b,c).
\end{aligned}
\]

对五个累计顶点 (p_0=e,p_1=g_1,p_2=g_1g_2,p_3=g_1g_2g_3,p_4=g_1g_2g_3g_4)，令

\[
 Z=T(p_1,p_2,p_3,p_4)-T(p_0,p_2,p_3,p_4)
 +T(p_0,p_1,p_3,p_4)-T(p_0,p_1,p_2,p_4)+T(p_0,p_1,p_2,p_3).
\]

逐面相消给出 (partial Z=0)。因为 (dV_{S^3}) 的周期是 (2\pi^2\)，

\[
 \int_Z dV_{S^3}\in2\pi^2\mathbb Z,qquad
 (\delta\omega_k^A)(g_1,g_2,g_3,g_4)=1.
\]

这个证明覆盖重复顶点、对跖点、(D=0) 和凸包包含原点的分支；随机五边形只能作为实现检查。

### Q8 精确核验

取

\[
 a=i\sigma_1=(0,1,0,0),\qquad b=i\sigma_2=(0,0,1,0),
\]

满足 (a^4=1,b^2=a^2,bab^{-1}=a^{-1})。程序 `scripts/ade_subgroups.py` 用有理四元数精确生成 Q8；不能把 `cos(\pi/2)` 的浮点近似当成精确输入，否则会把退化分支误判为一般位置。

`scripts/check_ade_review.py` 对 Q8 的全部

\[
 |Q_8|^4=4096
\]

组 ((g_1,g_2,g_3,g_4)) 检查 (deltaomega_1^A=1)。实际最大复数残差为

\[
 1.80\times10^{-86}.
\]

取 (C_4=\langle a\rangle) 的规范化 bar 三循环，按 A 的符号得到

\[
 P_4(\omega_1^A)=\prod_{r=0}^3\omega_1^A(a,a^r,a)=+i.
\]

这项数值不是 Q8 类阶为 4 的证明；它只核对符号和 (C_4) 限制。Q8 的完整类阶由第 4 节的有限子群定理识别。

## 3. (C_n) 限制与 (k=1) 非恰当性

取 (h=e^{2\pi i\sigma_1/n})，(0\le a<n)，并定义进位函数

\[
 N_n(a,b)=\left\lfloor\frac{a+b}{n}\right\rfloor,
 \qquad [a+b]_n=a+b-nN_n(a,b).
\]

标准 bar 三循环为

\[
 z_n=\sum_{r=0}^{n-1}[h\mid h^r\mid h].
\]

它的边界望远镜相消。A 约定的限制代表为

\[
 \boxed{
 \omega_{k,C_n}(h^a,h^b,h^c)
 =\exp\!\left(\frac{2\pi i k a}{n}N_n(b,c)\right).}
\]

直接配对得到

\[
 P_n(\omega_{k,C_n})=
 \prod_{r=0}^{n-1}\omega_{k,C_n}(h,h^r,h)
 =e^{2\pi i k/n}.
\]

若 (omega_1=\delta\beta) 在整个 SU(2) 上恰当，则限制到 (C_n) 仍为 (delta(\beta|_{C_n}))，其对 bar 三循环的乘积必须为 1；取 (n=4) 得 (P_4=i\ne1)。因此

\[
 \boxed{[\omega_1]\ne0;\quad \omega_1\text{ closed but not exact}.}
\]

这一步使用的是具体代表的子群配对，不是仅凭 (H^4(BSU(2),\mathbb Z)=\mathbb Z) 的抽象存在性论证。乘上任意二余边界不会改变 (P_n)。

## 4. (2D_n=\mathrm{Dic}_n)：精确群律和同一解析求值规则

固定

\[
 a=e^{i\pi\sigma_1/n},\qquad b=i\sigma_2,
\]

并以

\[
 x=(r,\varepsilon)\equiv a^r b^\varepsilon,
 \quad r\in\mathbb Z/(2n),\quad \varepsilon\in\{0,1\}
\]

标记元素。精确乘法是

\[
 (r,\varepsilon)(s,\delta)=
 \left(r+(-1)^\varepsilon s+n\varepsilon\delta\pmod{2n},
 \varepsilon+\delta\pmod2\right).
\]

对应四元数是

\[
 q_{r,0}=\left(\cos\frac{\pi r}{n},\sin\frac{\pi r}{n},0,0\right),
 \quad
 q_{r,1}=\left(0,0,\cos\frac{\pi r}{n},\sin\frac{\pi r}{n}\right).
\]

给定三个标签 (x_j=(r_j,\varepsilon_j))，用上述群律计算累计标签

\[
 p_0=(0,0),\quad p_1=x_1,\quad p_2=x_1x_2,\quad p_3=x_1x_2x_3,
\]

再把 (p_j) 映射到 (q_{r,\varepsilon})，形成 (Q,\Gamma,D)。这给出了覆盖所有三元组的明确求值规则：

1. (D\ne0)：使用任务 1 的八项 (operatorname{Li}_2) 公式；源 A 约定在指数前取负号。
2. (D=0)：用同一确定的 (E/F/T) 锥填充；每个面至多 6 项、候选至多 19 个，三维填充至多 24 个四面体、候选至多 73 个。每个非退化项仍由同一八项 (operatorname{Li}_2) 公式求值，合法降秩项的三形式积分为零。

因此这不是“把群元代入一个未定义体积函数”，而是一个有限、分支已固定、覆盖所有 (\mathrm{Dic}_n^3) 的算法。由于 (n) 任意，不能把所有三元组打印成表；标签规则和有限公式共同定义每个输入的相位。

`python -m scripts.check_ade_review` 的实际结果：

| 检查 | 结果 |
|---|---:|
| (mathrm{Dic}_n) 群律闭合，(n=2,\ldots,8) | 0 failures |
| (2T,2O,2I) 元素计数 | 24, 48, 120 |
| (2T,2O,2I) 精确代数乘法闭合 | 24, 48, 120 |
| Q8 全部四元组闭性 | (4096) 组，最大残差 (1.80\times10^{-86}) |
| Q8 (C_4) 配对 | (+i)（A 约定） |

## 5. 全部 ADE 的类阶和 (k) 周期

附件中给出的 ADE 分类为

\[
\begin{array}{c|c|c}
\text{类型}&\Gamma\subset SU(2)&|\Gamma|\\ \hline
A_{n-1}&C_n&n\\
D_{n+2}&\mathrm{Dic}_n&4n\\
E_6&2T&24\\
E_7&2O&48\\
E_8&2I&120
\end{array}
\]

周期群的标准结果是

\[
 H^3(\Gamma,U(1))\cong H^4(B\Gamma,\mathbb Z)\cong\mathbb Z_{|\Gamma|}
\]

对二元循环、二元二面体及二元多面体群成立；本次的 (2T,2O,2I) 精确 quaternion 集合和闭包已由脚本核验。Epa--Ganter, arXiv:1605.09192v1, Theorem 1.1 明确指出：String(3) 限制到任意有限 (\Gamma\subset S^3=SU(2)) 的类具有精确阶 (|\Gamma|)。因此本题的 CS 整类（A 约定只改变其逆元符号）满足

\[
 \operatorname{ord}[\omega_k|_\Gamma]
 =\frac{|\Gamma|}{\gcd(k,|\Gamma|)},
 \qquad
 [\omega_k|_\Gamma]=0
 \Longleftrightarrow |\Gamma|\mid k.
\]

对应周期为

\[
\begin{array}{c|c|c}
\Gamma&\text{类周期}&\text{平凡条件}\\ \hline
C_n&n&n\mid k\\
\mathrm{Dic}_n&4n&4n\mid k\\
2T&24&24\mid k\\
2O&48&48\mid k\\
2I&120&120\mid k
\end{array}
\]

这里的 ([\omega_k|_\Gamma]^{|\Gamma|}=0) 是上同调类的陈述；不表示未经二上链变换的几何代表逐点满足 (omega_k(g,h,l)^{|\Gamma|}=1)。也不表示任何具体物理体系的所有 gauging 条件已自动满足。

## 6. 这一轮的结论和未完成项

**已经核验：** 新笔记 A 的整体号与旧仓库相反，已给出精确转换；`src/su2_omega.py` 现在同时提供旧代表和 `omega_quaternions_source` 的 A 约定入口；全局闭性证明适用于 Q8 等异常分支；Q8 全部 4096 个四元组已用精确有理四元数核验；(C_n) 配对给出 (k=1) 非恰当性；Dic_n 精确群律和 (2T,2O,2I) 精确元素闭包已实现；全部 ADE 的类阶、周期和平凡条件由有限子群 String(3) 定理识别。

**尚未声称：** 尚未为 (mathrm{Dic}_n) 和 E 型群逐个打印一个短的初等三余循环表达式；当前交付是覆盖所有三元组的有限 Li₂/EFT 算法，以及由特征类定理识别的类。下一步若继续，应首先把 `source_phase` 的整体号转换正式并入主接口，再增加 Dic_n 的代表性非交换三元组数值表；不应把浮点三角值当作精确群关系。

## 7. 可复现命令

```bash
python -m scripts.check_ade_review
```

输出 JSON：`docs/review-0922/SU2_ADE_validation.json`。
