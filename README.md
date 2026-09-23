# SU(2) 对称缺陷三余循环

从 1+1 维玻色体系的 SU(2) 反常内流出发，给出一般群元的有限二重对数表达，以及覆盖退化输入的全局有限填充规则。

**阅读入口：[完整综合 PDF](paper/su2_complete_report.pdf) · [主研究笔记 PDF](paper/su2_cocycle.pdf) · [ADE 补充报告 PDF](paper/su2_ade_review.pdf) · [ADE LaTeX 源码](paper/su2_ade_review.tex)**

## 当前结果（2026-09-23 第二轮，分支 `polar-dual-exact`）

[第二轮完整推导](docs/review-0923/SU2_POLAR_DUAL_EXACT.md)：极对偶 + Schläfli 公式 + Coxeter 室，
给出**同一个四面体**的不含 Li₂ 的精确体积。本轮不换代表、不做陪集转移，
推导与证书中不用数值拟合或数值归约。

| 范围 | 第二轮结果 |
|---|---|
| 2O/2I 非退化类型 | 84/563 类全部有精确初等闭式 V=aπ²+(π/2)Σb_jΘ_j；2O 的 Θ=arctan(1/√2)，2I 的 Θ=arctan(1/2)、arctan(1/√5)、arctan√15 |
| 有理体积 | 2T 12/12、2O 60/84、2I 268/563；其余类型的无理性由线性无关性定理严格证明 |
| 群元直接相位 | 六个迹、取向迹和一个正负号室计数 N；无需按类型查表，只用角度证书表 |
| E 型退化输入 | 保持原 E/F/T 锥点，化为群四面体加月牙、悬挂（及球面型）初等修正；2T/2O/2I 全部 1852416 个有序输入逐一精确全检，失败 0 |
| D 表 | 在连续群 K=N(U(1)) 上形式证明闭性与归一化，对任意 n 成立；与几何代表同类，类阶 4n |
| 类阶 | 度数引理：几何代表在自由作用的有限群 G 上的类阶为 \|G\|（2T/2O/2I 为 24/48/120，Dic_n 为 4n） |

退化输入的相位依赖原锥点 a(1,t,t²,t³) 的坐标，只能写成与锥点有关的初等角，不能归入有限类型表。
需要专家重点审阅的是：Schläfli 公式、Humphreys §1.12/§1.14、Wigner 定理与 Serre 谱序列，以及度数引理。
综合 PDF 的第 7 节是证明与证书，第 8 节是全部类型的短表。本文是专家审阅稿，未经外部同行评议。

```python
from sympy import Rational, sqrt
from su2_polar_dual import polar_formula, polar_catalogue
from su2_polar_global import global_polar_formula
from su2_dic_exact import formal_closure_certificate, dic_pairing

h = Rational(1, 2)
r = polar_formula('2O', (0, 1, 0, 0), (1/sqrt(2), 0, 0, 1/sqrt(2)), (h, h, -h, h))
r['volume']            # -pi**2/16 + pi*atan(sqrt(2)/2)/2，即 (π/2)(β−π/8)，精确值
m = (-1, 0, 0, 0)
global_polar_formula('2T', m, m, m)['phase']   # -1：退化输入，保持原锥点
dic_pairing(7)         # Fraction(1, 28)：D 表与 S^3/Dic_7 基本类的配对 1/(4n)
```

输入必须是精确坐标（SymPy 有理数与根式）；模块拒绝浮点输入。

```bash
python -m unittest tests.test_polar_dual tests.test_polar_global tests.test_dic_exact -v   # 精确测试
PYTHONPATH=src python -m scripts.check_polar_dual                                       # 精确目录与证书 JSON
PYTHONPATH=src python -m scripts.check_polar_global_existence --group 2O --output results/x.json  # 逐输入全检
PYTHONPATH=src python -m scripts.render_exact_catalogue                                 # PDF 短表与证书表
PYTHONPATH=src python -m scripts.audit_polar_numeric                                    # 数值审计（不参与证明）
```

也可用 `make polarcheck`、`make polarreport`、`make polaraudit`。

## 第一轮路线（2026-09-23）

[第一轮推导](docs/review-0923/SU2_LI2_ADE_EXACT.md)：完整 Li₂ 解析式 → 纯符号恒等式 → 同一个四面体的直接几何。

| 范围 | 第一轮结果 |
|---|---|
| 一般 SU(2) | 三个群元的六个迹、累计取向迹、完整主支 Li₂/log 程序 |
| 纯符号化简 | 正交例、任意双圆弧族的完整证明；带适用域的共轭/Euler/倍角约简 |
| 2T | 12 种非退化无向类型，F4 几何证明及实际有理体积短表 |
| 2O/2I | 全部 84/563 种非退化类型的精确参数与解析 Li₂ 相位，无理体积反例 |
| E 型退化输入 | 保持原 E/F/T，至多 24 项的完整精确解析 Li₂ 表达 |

第一轮的逐类 CAS 化简只完成了一部分（659 类中 53/91/515），现已被第二轮的极对偶闭式完全取代；
其中已化简的 144 条记录经数值审计，与新闭式一致。

### 纯符号使用

Python 3.11+、SymPy 1.14.0：`python -m pip install -e '.[exact]'`。

```python
from sympy import pi
from su2_symbolic import group_formula
from su2_2t_geometry import phase_2t
from su2_exact_ade import global_formula, classify_group

triple = ((0, 1, 0, 0), (0, 0, 0, 1), (0, 1, 0, 0))
formula = group_formula(*triple)
assert formula['raw_volume'] == pi**2 / 8
geometry = phase_2t(*triple)  # 同一四面体，72 个 F4 室
types = classify_group('2I')  # 563 行，返回独立副本
center = (-1, 0, 0, 0)
full = global_formula('2I', center, center, center)  # 完整精确表达，不强制 CAS 化简
```

新边序为 `(01,02,03,23,13,12)`，主式返回体积模 `2*pi**2`。
入口拒绝 Float、非单位群元、非整数 level 和非法群成员。
单四面体入口拒绝退化；全局入口处理退化。
只需完整通用式时可用 `group_formula(..., reduce=False)`；默认符号约简可能较耗时。
`expand=False` 仅作链审计，有未展开的满秩项时明确返回 `phase=None`。

### 完整参数表和精确复现

[精确 JSON](results/SU2_symbolic_ADE_exact.json) 保存全部 659 类型的代表、Gram 键、重数、
六个余弦、行列式和完整公共 `formula_template`；每行 `formula_inputs` 可直接代入。
PDF 内也印有全部 84/563 个 Gram 参数键及解码规则，专家无需运行程序才能读到完整表。
另含 2T 体积表、两组无理性证书和逐类型符号尝试。
CAS 超时仅表示该次尝试未完成；缓存与耗时会影响尝试状态，不影响精确参数表。

```bash
# 新路线测试：纯精确，不作浮点拟合
python -m unittest tests.test_su2_symbolic tests.test_2t_geometry tests.test_exact_ade tests.test_exact_certificates tests.test_symbolic_cli -v

# 完整参数/公式表，无需执行耗时的逐类尝试
python -m scripts.check_symbolic_2t --symbolic-trials none --output build/exact_catalogue.json

# 每个类型都尝试纯符号化简，并保存结果或明确的超时状态
python -m scripts.check_symbolic_2t --symbolic-trials all --symbolic-budget-seconds 4 --jobs 4 --output results/SU2_symbolic_ADE_exact.json

# 任意类型的完整展开式，不作浮点求值
python -m scripts.check_symbolic_2t --group 2I --type-index 0 --output build/2I_type0.json
python scripts/build_complete_report.py --engine tectonic
```

也可用 `make symboliccheck`、`make symbolicreport`、`make completepdf`。
完整兼容性回归需安装 `.[validation]`；其中旧 mpmath 测试不作为本轮数学证据。
旧 `make verify` 会执行历史数值检查，不是新纯精确路线的必要步骤。

第二轮代码：[极对偶体积](src/su2_polar_dual.py)、[退化输入初等相位](src/su2_polar_global.py)、[D 表任意 n 证明](src/su2_dic_exact.py)。
第一轮代码：[Li₂ 主程序](src/su2_symbolic.py)、[2T 几何](src/su2_2t_geometry.py)、
[三群分类及全局相位](src/su2_exact_ade.py)、[本轮 PDF 章节](paper/su2_symbolic_ade.tex)。
完整参数与体积短表由 `PYTHONPATH=src python -m scripts.render_exact_catalogue` 纯精确枚举生成。

## 历史材料（以下截至 09-22）

以下说明保留旧接口及当时的完成度，当前范围以上面的 09-23 说明为准。
[旧计数代表研究](docs/review-0922/SU2_E_ALGEBRAIC_RESEARCH.md) 给出同类但逐点不同的代表及比较余链，
**不作为本轮用户要求的解法**。
[历史人类 review 答复](docs/review-0922/SU2_HUMAN_REVIEW_RESPONSE.md) 保留当时的边长公式审计。

采用反厄米联络、基本表示普通迹，固定

\[
I_4=-\frac{k}{8\pi^2}\operatorname{Tr}(F\wedge F),\qquad k\in\mathbb Z.
\]

一般位置的相位通过八个主支二重对数求值；全局代表通过预先相容的边、面和确定锥点定义，最多包含 24 个球面四面体。文稿给出五边形、归一化、物理二级类识别及约定改变的证明。它是选定结点约定下的一个代表。

## 仓库结构

```text
paper/          自包含 LaTeX 文稿及编译后的 PDF
src/            主算法 su2_omega.py
scripts/        独立积分、级数、棱长复算和文稿构建
tests/          单元、解析校准及链结构测试
results/        固定样本、参数、精度和实际误差的 JSON 报告
docs/notes/     详细推导及复核笔记
docs/history/   原始交接、执行计划和历史审计记录
references/     用户提供的双语论文参考件及来源说明
```

## 安装与计算

完整复现环境为 Python 3.13、mpmath 1.3.0、NumPy 2.4.4、SymPy 1.14.0；主算法只需要 mpmath，E 型精确代数检查还需要 SymPy。建议在独立虚拟环境中安装：

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[validation]'
```

```python
import numpy as np
from su2_omega import omega_su2, omega_quaternions

center = -np.eye(2, dtype=complex)
value = omega_su2(center, center, center, k=1, dps=50)  # ≈ -1

# 四元数按正范数单位化。整数／Fraction 输入可精确保留群关系。
exact_input_value = omega_quaternions(
    (-1, 0, 0, 0), (-1, 0, 0, 0), (-1, 0, 0, 0), k=1
)
```

直接使用边长路线（A 负号约定），以及 E 型完整计算记录：

```python
from su2_omega import edge_lengths_quaternions, omega_quaternions_edge
from scripts.ade_phase import phase_details  # 从仓库根目录运行
from sympy import Rational

triple = ((0, 1, 0, 0), (0, 0, 0, 1),
          (Rational(1, 2), Rational(1, 2), -Rational(1, 2), -Rational(1, 2)))
lengths = edge_lengths_quaternions(*triple)  # 次序 23,13,03,01,02,12
phase = omega_quaternions_edge(*triple)     # exp(-i*pi/32)
details = phase_details(*triple)           # 精确输入、分支、六边长、体积和相位
```

新增 E 型代数代表（与几何代表同类，但逐点不同）：

```python
from scripts.ade_count import count_details, algebraic_phase, gauge_details

exact = count_details('2T', *triple, k=1)
assert exact['count'] == 1
assert exact['exponent_mod_order'] == 23  # exp(2*pi*i*23/24) = exp(-i*pi/12)
phase_alg = algebraic_phase('2T', *triple, dps=60)
beta = gauge_details('2T', triple[0], triple[1], dps=45)
# 对全部三元组：phase_alg = phase_geom * delta(beta)。完整证明见上面的研究稿。
```

`count_details` 保留未约化的有向整数计数和每个贡献者的精确群元；相位指数才按群阶取模。入口拒绝浮点坐标、非群成员和非整数 level。`gauge_details` 的体积比较使用高精度数值计算，返回全部平均项；它不参与整数计数的分支判断。

若需要乘积恰为单位元或中心元，使用 `quaternion_multiply` 在精确四元数表示中完成乘法。矩阵入口检查尺寸和 SU(2) 数值残差，并将允许的舍入残差投影到明确的归一化四元数。它不会恢复浮点输入背后未被表示的精确代数关系。

## 重现检查

在仓库根目录运行：

```bash
make verify
```

等价的逐项命令为：

```bash
python -m unittest discover -s tests -v
python -m scripts.su2_audit_checks
python -m scripts.validate_su2
python -m scripts.su2_edge_crosscheck
python -m scripts.check_ade_phase
python -m scripts.check_e_algebraic
```

主验证先生成 `results/SU2_validation.json`，棱长路线随后读取其中固定的样本并生成 `results/SU2_edge_crosscheck.json`。报告路径由脚本位置确定。`requirements.txt` 固定完整计算环境中的两个数值包版本。

| 检查 | 已执行结果 |
|---|---|
| 单元测试、精确链边界、异常输入 | 39 项通过；含 E 型字段、18 个分支样例的边界/等变检查和代表转换 |
| 24 个一般输入与原三重积分 | 最大体积差约 `1.85e-13`，正负取向各 12 个 |
| 棱长闭式独立复算 | 原固定样本扩展到全部 24 个，报告同时比较生产边长、独立边长及二面角路线 |
| E 型逐项检查 | 每群 3 组五边形；保存全部五相位及其计算项，另有手算例和退化校准 |
| 新 E 型代数代表 | 每群 16 组整数五边形，模群阶的余数均严格为 0；C4 配对为 +i |
| 显式 β 转换 | 每群一个三元组、四个 β，共 768 个平均项；45 位输出的最大残差约 6.12e-46 |
| 全局校准 | 中心值 `(-1)^k`、C4/C8 不变量、五边形、归一化均通过 |

这些误差为实际交叉检查结果，未提供任意输入下的区间算术误差证书。数学证明和数值检查各自承担不同的验证作用。

## 构建 PDF

使用 Tectonic 0.15.0 和系统字体 **Droid Sans Fallback、DejaVu Serif**。Tectonic 首次运行可能下载 TeX 依赖；可用 `fc-match 'Droid Sans Fallback'` 检查中文字体。

```bash
make pdf

# ADE 补充报告
make adepdf

# 编译单一源入口的完整专家审阅稿（不是拼接两个 PDF）
make completepdf
```

原推导的 LaTeX 源文件是 `paper/su2_cocycle.tex`；完整报告入口为 `paper/su2_complete_report.tex`，复用原推导与 ADE/E 型章节，统一编译。构建脚本检查缺字、未解析引用和溢出的排版盒子，通过后更新对应 PDF。`make all` 顺序运行全部检查和 PDF 构建。

## 推导与来源

- [完整数学推导](docs/notes/SU2_derivation.md)
- [全局链证明](docs/notes/SU2_global_review.md)
- [Chern–Simons 类及符号](docs/notes/SU2_physics_review.md)
- [完成度二次复核及独立棱长公式](docs/notes/SU2_completion_doublecheck.md)
- Jia 等：[arXiv:2510.14722v2](https://arxiv.org/abs/2510.14722v2)
- Murakami：[arXiv:1011.2584v4](https://arxiv.org/abs/1011.2584v4)
- Baez–Lauda：[arXiv:math/0307200v3](https://arxiv.org/abs/math/0307200v3)

历史审计记录保留当时的测试数量、缺陷及修复过程；当前运行入口、结果和文稿以上述文件为准。代理辅助审查记录属于研究过程材料，不表示外部同行评审。
