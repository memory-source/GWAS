## Round 15 (run_20261007_66850dac) — Score: 8.5/10 — regression check

### Raw Reviewer Response (verbatim, summary)
第 15 轮回归检查：W13-1..W13-5 全部通过产物验证；规模匹配重跑 s48_scalematched.json 已验证（3pop +0.161 Fst 0.0997/n=60；5pop +0.111 Fst 0.0885/n=100；单人群 +0.091/+0.076；控制 0.518–0.562）。评分 8.5（恢复）。结论：是 — 准备好提交。
剩余润色项(非阻碍): chr2 3pop 异常解释 "consistent with its small size" 不实(chr2 3pop n=450 与 chr1 相同; 正确解释 proximity to crossing); "n_study=20–36" 应为 28–36; s49 控制范围 "0.53–0.60" 应为 0.55–0.60; 引言 "below the real-data level" 早于共同框架。

### Memory Update
- 周期总结：R11 7.5→8.5 → R12 8.5 → R13 定律轮 8.0 (发现 W13-1..5) → R15 8.5 是(规模匹配重跑暴露混淆; 重构为"共同受规模和结构控制").
- 未来审查方法：diff 快照; 按 JSON 键枚举引用范围全部点; "monotone" 检查相邻梯级; 结构轴检查共变规模; ≈1.0 AUROC 算 k/n.