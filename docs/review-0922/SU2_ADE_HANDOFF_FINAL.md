# SU(2) ADE 任务最终交接审计

日期：2026-09-22。本文逐项对应 `SU2_anomaly_ADE_agent_prompt_zh(1).md`，并记录本轮发现的过早完成声明、修复、实际命令和剩余边界。提交前必须读本文；它比早期 `SU2_ADE_progress.md` 更严格。

## 结论先行

上一轮有一项不够严谨：E 型的相位入口虽然能跑代表性样本，但没有在二次域乘法后统一约化、没有可靠处理代数正负号，且锥点归一化引入了不必要的扩域。因此当时“E 型已完成”的措辞过强。本轮已经修复这三个问题，并重跑 E 型代表性五边形检查。

现在的状态是：

- 任务 1 的完整一般位置八项 Li₂ 公式（不是只给链接）在 `paper/su2_cocycle.pdf` 第 4 节，且已逐项复制到 `paper/su2_ade_review.pdf` 的开头；A 源笔记负号由 `omega_quaternions_source` 明确实现。
- 任务 2 的 closed 证明、(C_4) 非恰当性和 Q8 全部 4096 组精确检查已完成。
- 任务 3 的 (C_n)、(mathrm{Dic}_n)、(2T,2O,2I) 均有确定的有限求值入口；Dic 表对 (n=2,3,4,5) 全部四元组做了有理数闭性检查；E 型入口已使用精确二次域分支判定并通过代表性五边形检查。
- ADE 类阶、(k) 周期和平凡条件由有限子群 String(3) 定理识别；这一步不依赖把有限群的所有三元组打印出来。

参考附件没有被提交到 Git；它们仍在本地 `docs/review-0922/`。当前远端最新提交见仓库 `main`。

## 任务 1：最终公式确切位置和内容

一般位置 (D\neq0) 的最终公式在：

- `paper/su2_cocycle.tex` / `paper/su2_cocycle.pdf` 第 4 节；
- `paper/su2_ade_review.tex` / `paper/su2_ade_review.pdf` 第 2 页附近的“任务 1 的最终可代入公式”；
- `docs/notes/SU2_derivation.md` §1。

输入是
\[
p_0=1,\quad p_1=g_1,\quad p_2=g_1g_2,\quad p_3=g_1g_2g_3,
\qquad Q=(p_0,p_1,p_2,p_3),\quad D=\det Q,\quad \Gamma=Q^TQ.
\]
令 (H=\Gamma^{-1})，六个面对按
\[
(01;23),(02;13),(12;03),(23;01),(13;02),(03;12)
\]
编号，其中前一对是相对面，后一对是实际棱。二面角为
\[
\theta_j=\arccos\left(-\frac{H_{u_jv_j}}
{\sqrt{H_{u_ju_j}H_{v_jv_j}}}\right),qquad a_j=e^{i	heta_j}.
\]
辅助量为
\[
\begin{aligned}
r_0={}&a_1a_4+a_2a_5+a_3a_6+a_1a_2a_6+a_1a_3a_5+a_2a_3a_4\\
&+a_4a_5a_6+a_1a_2a_3a_4a_5a_6,\\
r_1={}&4\sum_{j=1}^3\sin\theta_j\sin\theta_{j+3},qquad r_2=\bar r_0,\\
z={}&-\frac{2r_0}{r_1+4\sqrt{\det\mathsf G}},qquad |z|<1.
\end{aligned}
\]
八项二重对数组合是
\[
\begin{aligned}
\mathcal L=\frac12\bigg[&\Li_2(z)+\Li_2\left(\frac z{a_1a_2a_4a_5}\right)
+\Li_2\left(\frac z{a_1a_3a_4a_6}\right)
+\Li_2\left(\frac z{a_2a_3a_5a_6}\right)\\
&-\Li_2\left(\frac{-z}{a_1a_2a_3}\right)
-\Li_2\left(\frac{-z}{a_1a_5a_6}\right)
-\Li_2\left(\frac{-z}{a_2a_4a_6}\right)
-\Li_2\left(\frac{-z}{a_3a_4a_5}\right)
-\sum_{j=1}^3\theta_j\theta_{j+3}\bigg],\\
W={}&-\Re\mathcal L+\pi\left(\Arg(-r_2)+\frac12\sum_{j=1}^6\theta_j\right)-\frac32\pi^2.
\end{aligned}
\]
Murakami 定理给 (V=W\bmod2\pi^2)。A 源笔记的网络比值约定使用
\[
\boxed{\omega_{k,A}=\exp\left(-\frac{ik\operatorname{sgn}(D)}\pi W\right)}.
\]
平方根取正实根，(Argin(-pi,pi])，(Li_2) 取主支。(D=0) 由 E/F/T 有限填充处理；每项仍调用同一八项公式。

## 任务 2：closed 与 non-exact

E/F/T 链满足
\[
\partial T(a,b,c,d)=F(b,c,d)-F(a,c,d)+F(a,b,d)-F(a,b,c).
\]
五个累计顶点组成的五项三链 (Z) 逐面相消，所以 (partial Z=0)。由于 (dV_{S^3}) 的周期是 (2\pi^2)，
\[
\delta\omega_{k,A}=1
\]
对所有异常输入成立。该证明是链级证明，随机数值只作实现检查。

在 (C_n=\langle h\rangle)、(h=e^{2\pi i\sigma_1/n}) 上，进位代表为
\[
\omega_{k,C_n}(h^a,h^b,h^c)=exp\left(\frac{2\pi i k a}{n}\left\lfloor\frac{b+c}{n}\right\rfloor\right).
\]
规范化 bar 三循环 (z_n=\sum_r[h|h^r|h]) 给出
\[
P_n=\prod_r\omega_{k,C_n}(h,h^r,h)=e^{2\pi ik/n}.
\]
取 (n=4,k=1)，(P_4=i\ne1)，而余边界在三循环上的配对必为 1。因此 (omega_1) 非恰当。

Q8 的精确检查使用有理四元数，不使用浮点 (cos(pi/2))。全部 4096 个四元组的最大闭性残差为 (1.80\times10^{-86})。

## 任务 3：全部 ADE

### (C_n)

相位和周期如上，类阶为 (n/\gcd(k,n))。

### (mathrm{Dic}_n=2D_n)

精确标签为 (a^rb^\epsilon)，(a^{2n}=1,b^2=a^n,bab^{-1}=a^{-1})，乘法为
\[
(r,\epsilon)(s,\delta)=\bigl(r+(-1)^\epsilon s+n\epsilon\delta\pmod{2n},\epsilon+\delta\pmod2\bigr).
\]
另有源半球代表的八分支表，见 `scripts/dic_formula.py`。其有理数闭性在 (n=2,3,4,5) 的全部四元组上通过。任意 n 的证明使用统一边面链；表格中的余数满足
\[
[x]_m=x-m\left\lfloor\frac{x}{m}
ight
floor,
\qquad N(a,b)+N([a+b]_m,c)=N(b,c)+N(a,[b+c]_m),
\]
其中 m=2n。将八种 s_1s_2s_3 分支代入五项余边界，逐项用这个进位恒等式和 [x]_m 的定义约去，剩余项都是 2Z 中的整倍数。因此 exp(-i*pi*k*delta v_n)=1 对任意 n,k 成立；有限穷举只承担实现回归，不承担普遍证明。

八分支表与 E/F/T 几何代表可能使用不同的边链。严格比较不是一句“换面”：取
\[
\partial A(g)=\gamma_A(g)-\gamma_{EFT}(g),
\]
再取
\[
\partial B(g,h)=S_A(g,h)-S_{EFT}(g,h)-\bigl(gA(h)-A(gh)+A(g)\bigr).
\]
由于 H_1(S^3)=H_2(S^3)=0，这些比较链存在；于是
\[
\partial(C_A-C_{EFT}-\delta_GB)=0,
\qquad \omega_A=\omega_{EFT}\,\delta\beta.
\]
这才证明两套代表同类，而非逐点相等。

### (2T,2O,2I)

精确列表、二次域约化和 E/F/T 相位入口在：

- `scripts/ade_subgroups.py`：精确元素集合和闭包；
- `scripts/ade_phase.py`：精确分支、未归一化 moment-ray 锥点、Murakami 求值；
- `scripts/check_ade_phase.py`：代表性非交换五边形检查。

关键修复是：

1. 每次四元数乘法后按 (s_2^2=2) 或 (s_5^2=5) 约化；
2. 用正根代数符号判断凸包正依赖，不能使用形式符号的 `is_nonnegative`；
3. 锥点使用正射线 ((1,t,t^2,t^3))，不引入 (sqrt{85}) 扩域；
4. 只有分支确定后，才将代数数转换为 mpmath 求 Li₂ 体积。

这三个问题是本轮审计发现并修复的，不是仅仅排版调整。修复后 (2T,2O,2I) 各两组非交换样本的五边形最大残差分别为
\[
9.56\times10^{-67},\quad1.41\times10^{-66},\quad2.33\times10^{-66}.
\]

先由 SU(2) 上的 CS 微分特征自然性识别
\[
[\omega_{A,k}|_\Gamma]=\varepsilon k[\mathrm{String}(3)|_\Gamma],
\qquad \varepsilon=\pm1,
\]
其中 epsilon 只记录 A 与正体积生成元的整体取向。再用 Epa–Ganter Theorem 1.1 说明 String(3) 限制到有限 Gamma\subset SU(2) 的类精确阶为 |Gamma|。因此
\[
H^3(\Gamma,U(1))\cong H^4(B\Gamma,\mathbb Z)\cong\mathbb Z_{|\Gamma|},
\qquad
\operatorname{ord}[\omega_k|_\Gamma]=\frac{|\Gamma|}{\gcd(k,|\Gamma|)}.
\]
ADE 周期和平凡条件是：
\[
\begin{array}{c|c|c}
\Gamma&\text{周期}&\text{平凡条件}\\\hline
C_n&n&n\mid k\\
\mathrm{Dic}_n&4n&4n\mid k\\
2T&24&24\mid k\\
2O&48&48\mid k\\
2I&120&120\mid k
\end{array}
\]

## 实际文件和命令

```bash
python -m unittest discover -s tests -v
python -m scripts.check_ade_review
python -m scripts.check_dic_table
python -m scripts.check_ade_phase
make adepdf
```

综合 PDF 是 `paper/su2_cocycle.pdf`；本次 ADE 专题 PDF 是 `paper/su2_ade_review.pdf`。源码、JSON 报告和验证脚本已在远端 `main` 分支。用户提供的 PDF/JPG 参考附件仍是本地未跟踪文件，没有进入提交。

## 最终边界

现在“全部三元组可求值”指一个确定的、精确分支选择后调用有限 Li₂ 体积的程序入口，而不是已经生成所有 (2I^3) 项的静态表。类阶和平凡条件是定理级识别；数值残差只核查实现。两者没有混为一谈。
