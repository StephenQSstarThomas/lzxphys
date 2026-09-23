# 计算报告

- `SU2_validation.json`：`python -m scripts.validate_su2` 生成。固定随机种子 20260913，记录 24 个一般输入的精确有理四元数原始整数坐标、闭式体积、逐次提高阶数的独立三重积分、误差和各项检查。额外记录共轭、level 相加、五边形、有限子群、中心、近退化精度比较，以及两种性质不同的分支行为。
- `SU2_edge_crosscheck.json`：`python -m scripts.su2_edge_crosscheck` 生成。覆盖前一报告的全部 24 个样本，比较生产边长路线、独立边长实现和二面角路线；分别保存两种差值。
- `../docs/review-0922/SU2_ADE_phase_validation.json`：`python -m scripts.check_ade_phase` 生成。每个 E 型群三组五边形，保存全部五项的精确输入、累计点、有限填充、边长、体积和相位，另有独立手算例及退化校准。

这些报告由源码生成并纳入版本控制，便于核对文稿中的数字。运行 `make verify` 会按正确依赖顺序重建它们。报告的 `passed` 是已列门控的合取，不是数值误差的严格上界证明。

`formula_dps` 等字段是请求的十进制精度；主实现还会根据几何消去和整数 level 的位数增加工作精度。独立三重积分使用 NumPy float64；两个收敛阶数之间的差和最终与闭式的差都保留在报告中。
