# 人类 review 修正计划

> **For agentic workers:** Use superpowers:executing-plans for task-by-task implementation; user has already authorized corrections and advancement.

**Goal:** 交付三个群元直接给出六边长的体积公式，以及可逐项复算的 E 型相位公式。

**Architecture:** 保留现有全局相容链选择，提升边长闭式为公开求值路线；E 型在精确代数分支确定后使用该路线。文档先给完整公式，再给源码映射和可审计数据。

**Tech Stack:** Python、mpmath、SymPy、unittest、LaTeX/Tectonic。

**Spec:** 本轮用户 review：不能重复二面角公式冒充边长答案；不能以 E 型数值检查或程序入口代替群元素公式。

## Global Constraints

- 源 A 约定为 exp(-ikV/pi)，群乘法采用 +i sigma 的负叉积。
- 六边长只能决定无向体积；群元累计顶点的行列式决定取向。
- 退化输入必须保留相容边面链，不能一概返回零。
- 当前 checkout 没有 09-21 命名材料；已有 review 在 docs/review-0922，缺失的原附件不作为已重新核验的证据。
- 直接修正当前工作区；不提交或推送。测试和最终独立审查后交付可检查 diff。

## Review Focus

边编号/取向；精确退化与对跖输入；非整数 level；E 型根式与实际群元素；报告是否真正记录五项而非只记录一个残差。

### Task 1: 群元到边长求值

- [x] 在 tests/test_review_formulas.py 写正交四面体 V=pi²/8、镜像取向、直接边长顺序、中心退化的测试；先运行观察缺失接口失败。
- [x] src/su2_omega.py 增加 edge_lengths_quaternions、oriented_volume_from_edges、omega_quaternions_edge，完整展开定理 1.2 的固定 z 对数项；scripts/su2_edge_crosscheck.py 保持原独立实现作比较。
- [x] 运行新测试和原 unittest 全套，比较固定样本的边长与二面角两条路线。

### Task 2: E 型可审计相位与实际缺陷

- [x] 测试 E 型非整数 k 拒绝、三个手算群元例、D 型 a^r b 坐标乘法；先观察错误。
- [x] scripts/ade_phase.py 添加可审计求值数据（精确输入、累计点、分支、每项六边长、符号、体积、相位），默认使用边长路线，保留 angle 路线交叉校验。
- [x] 修复 scripts/ade_subgroups.py 的 D 型符号，pyproject.toml 声明 SymPy 验证依赖。
- [x] scripts/check_ade_phase.py 保存每组四元数、五个相位、左右乘积和残差；增加各群非共享元素例及退化例，明确覆盖范围。

### Task 3: 数学交付与审计结论

- [x] 新建 docs/review-0922/SU2_HUMAN_REVIEW_RESPONSE.md，包含完整边长公式、E6/E7/E8 元素集合与乘法、有限分支相位公式、三组手算校准及发现的问题。
- [x] 用同一 LaTeX 公式章节更新 ADE PDF 与综合 PDF；修正链边界漏参数和历史完成声明入口。
- [x] 更新 README 和验证入口，运行 unittest、边长交叉检查、E 型逐项检查及 PDF 构建。
- [x] 按 requesting-code-review 技能执行一次独立只读审查，处理重要问题，记录实际结果及未覆盖边界。

## 执行记录

初始审计：二面角公式仍在 ADE 主交付；边长已有独立附录实现但未把六迹显式连接到主答案；E 检查仅每群两组。另发现 dic_quaternion 的末坐标号错误、E 入口 int(k) 截断、A 约定倒数在默认精度下计算、验证依赖遗漏 SymPy、ADE 文稿边界式 F(a,b,c) 漏 c。

最终验证：26 项 unittest 通过；原 24 个独立积分全部收敛，最大差 1.84742e-13；24 个独立边长/二面角比较最大差 2.118e-83，生产边长路线也全部通过。E 型各三组五边形，最大差分别 1.143e-65、9.108e-66、8.464e-66，三组手算例、中心与归一化校准均通过。基准审计通过。

补充验证：check_ade_review 已完成三组 E 元素精确闭包，以及 Q8 全部 4096 个四元组（最大差 1.138e-86、C4 配对 +i）；报告已更新。

独立只读审查没有发现重要问题；另外检验 24 个随机四面体、1e-20/1e-60/1e-100 的微小四面体以及 k=10^100+1 的三个 E 型例子。审查指出的边长记号和分母分组已修正。其后补的 normalized_numeric 正根转换和 D 型非整数 level 拒绝均先复现再修正，由 26 项全套测试覆盖。

PDF：Tectonic 构建成功，未出现 overfull、缺字或未解析引用；查看了边长与 E 型公式的实际页面。环境缺少 make，执行了其等价命令。综合 PDF 将修订 ADE 答复置前。未提交、推送。

边界：E 型交付通用有限几何公式，专属短表尚未导出；旧任意 n 的 D 表符号化证明未补齐；原源 A 附件不可重新核验。当前工作区直接修正符合用户授权，不另外引入分支/提交操作。
