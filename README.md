# SU(2) 对称缺陷三余循环

从 1+1 维玻色体系的 SU(2) 反常内流出发，给出一般群元的有限二重对数表达，以及覆盖退化输入的全局有限填充规则。

**阅读入口：[完整综合 PDF](paper/su2_complete_report.pdf) · [主研究笔记 PDF](paper/su2_cocycle.pdf) · [ADE 补充报告 PDF](paper/su2_ade_review.pdf) · [ADE LaTeX 源码](paper/su2_ade_review.tex)**

**本轮 E 型推进：[精确代数公式、代表转换与验收记录](docs/review-0922/SU2_E_ALGEBRAIC_RESEARCH.md)。** 保留原几何公式，新增只用群元、行列式和精确正负号的轨道计数代表，并给出显式二余链 β，证明新旧代表共边界等价。综合 PDF 使用单一标题、目录、连续编号和统一的相位约定；保留完整 SU(2)/CS 推导。它是专家审阅稿，不表示已经获得外部同行评审。

[前一轮人类 review 答复](docs/review-0922/SU2_HUMAN_REVIEW_RESPONSE.md) 保留从三个群元的六个迹到边长体积的完整公式。E 型现已给出有界精确代数公式，但未化为类似 D 型八行表的专属最短表；任意 n 的 D 表完整证明仍是单列遗留项。

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
