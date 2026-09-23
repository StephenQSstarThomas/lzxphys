# 计算报告

## 第二轮（极对偶，精确）

- `SU2_polar_dual_exact.json`：`PYTHONPATH=src python -m scripts.check_polar_dual` 生成。
  659 类型的室数 N、有理部分 a、基底系数 b，2O/2I 的角度表与精确证书，D 表的形式闭性、
  归一化残差和透镜及 Dic_n 基本闭链配对。全部是精确运算，没有浮点。
- `SU2_polar_global_existence.json`：`scripts.check_polar_global_existence` 对全部有序三元组的精确全检
  （2T/2O/2I 合并），记录分支、修正模式和失败数。
- `SU2_polar_numeric_audit.json`：`scripts.audit_polar_numeric` 生成的**数值审计**，
  与历史 mpmath 实现和精确 Li₂ 链交叉核对，不参与任何证明。

## 历史报告

- `SU2_validation.json`：`python -m scripts.validate_su2` 生成。固定随机种子 20260913，记录 24 个一般输入的精确有理四元数原始整数坐标、闭式体积、逐次提高阶数的独立三重积分、误差和各项检查。额外记录共轭、level 相加、五边形、有限子群、中心、近退化精度比较，以及两种性质不同的分支行为。
- `SU2_edge_crosscheck.json`：`python -m scripts.su2_edge_crosscheck` 生成。覆盖前一报告的全部 24 个样本，比较生产边长路线、独立边长实现和二面角路线；分别保存两种差值。
- `../docs/review-0922/SU2_ADE_phase_validation.json`：`python -m scripts.check_ade_phase` 生成。每个 E 型群三组五边形，保存全部五项的精确输入、累计点、有限填充、边长、体积和相位，另有独立手算例及退化校准。

这些报告由源码生成并纳入版本控制，便于核对文稿中的数字。运行 `make verify` 会按正确依赖顺序重建它们。报告的 `passed` 是已列门控的合取，不是数值误差的严格上界证明。

`formula_dps` 等字段是请求的十进制精度；主实现还会根据几何消去和整数 level 的位数增加工作精度。独立三重积分使用 NumPy float64；两个收敛阶数之间的差和最终与闭式的差都保留在报告中。
