## Round 2 (run_20261006_f68b0cd6) — Score: 8.0/10

### Raw Reviewer Response (verbatim)
第二轮审查已完成。结论已写入 `refine-logs/pi_review_verdict_round11_loop2.md`.

## W1–W8 评估（已与制品比对验证）
- W1 倒置 A* vs 2AUC−1: 已修复. 论文现为 "A* 与二项正态值 2AUC−1 不同（对于泄露通道低于该值，对于随机概率控制高于该值）" — 与 JSON 完全匹配（泄露差距为 −0.056…−0.153；控制为 +0.029/+0.043）.
- W2 ε_lb 公式/剪裁/n_neg: 已修复. 代码接收 n_pos,n_neg，在 1/n_pos 处对 FNR 进行下限处理；论文陈述了双分支最大值、1/900 分辨率以及排序信息警告；JSON 重新运行 pc 4.997→≈5.0, af_1dp 1.057→1.1；评估器合法.
- W3 未来工作冲突: 已修复. §replication 现在指向 §scale.
- W4 交叉引用: 已修复. 局限性 → sec:scale.
- W5 饱和扫描过度销售: 已修复. "使通道在 ρ=1 时保持饱和".
- W6 "full marker density": 部分修复. 括号披露了 "uniformly strided from the full-density parse"，但主句仍写 "at full marker density"，自相矛盾.
- W7 设计文档 / CV: 大部分修复. ROUND11_DESIGN.md 记录了交付的 6 单元网格；留下两处新的文档残留 (R2, R3).
- W8 死导入 / docstring: 已修复.

回归检查：六个预期修改块；s39/s40/s41/s42/s37 JSON 未受影响；s38 A*/AUC/gap 与第一轮相同（仅 eps_lb 变化），符合预测；编译成功 15 页.

## 剩余缺陷（按排序，含最小修复方案）
1. R1 (次要, 论文) — "versus 0.53/4.2 for full-precision AF" 引用了 af_2dp 的值 (0.5318/4.190)；af_full 为 0.5456/4.120 → 0.55/4.1.
2. R2 (次要, 文档) — B5 "donor CoV 186/141=1.32 vs cohort 50/159=0.32" 归因错误：186/141 是 nsp150 cohort 的 CV（donor 48/128=0.38），50/159 是 nsp10000 cohort 的 CV.
3. R3 (次要, 文档) — EXPERIMENT_TRACKER.md S38 行仍带旧方向速记 "A*−(2AUC−1)>0".
4. R4 (次要, 论文) — s41 主句仍写 "at full marker density".
5. R5 (可选) — 四个泄露通道的 eps_lb_0.01_valid 为 False（零 FN 点处严格检查不兼容任何有限 ε）；加一个说明性条款.

## 评分：8.0/10 → 在完成 R1–R4 后达到 8.5
## 结论：差不多 — 只差一行编辑的通过.

### Memory Update
- W1–W5, W8 完全修复；W6 部分修复；W7 大部分修复. 评分 8.0 → R1–R4 后 8.5.
- R1–R4 修复列表 (见上). R5 可选：披露 eps_lb_0.01_valid=False 标志.
- 验证稳定：新的 ε 评估器合法且与 g_{ε,δ} 一致；s38 A*/AUC 两次运行相同；其他 JSON 未受影响；论文 diff = 六个预期块；15 页编译成功.

---

## Round 3 (run_20261006_f68b0cd6) — Score: 8.5/10