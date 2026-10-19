# GWAS 项目独立审稿：2026-10-15

**按 NeurIPS / ICML / ICLR 主会研究论文评估，当前版本我会倾向拒稿。** 核心 PCA 成员泄漏现象有证据，主要模拟数值能从归档重算；决定性缺口是新贡献的力度、真实发布场景与攻击者能力的对应，以及在公平强基线之上的新增泄漏。继续增加同类网格和随机种子，难以解决这些问题。

这是一份对当前稿件和项目证据的独立评估，不是会议正式决定。目标会场尚未指定，暂采用上述机器学习主会的贡献、可靠性和可迁移性要求。安全隐私或计算生物学会场会有不同侧重，但同样需要修复下述发布对象和比较问题。报告采用三个审稿侧重点，三者共享同一组可核对事实，不虚构审稿人身份。

## 评审范围与验证边界

当前论文是 [paper/main.tex](/data/lvvyongji1/GWAS/paper/main.tex:13) / [paper/main.pdf](/data/lvvyongji1/GWAS/paper/main.pdf)，题为 *The Diagnostics Release Is the Leak: Controlled Measurement of Revealed Statistics in Multi-Site GWAS*。PDF 为 22 页，创建时间为 2026-10-14 11:55:29 +08:00。TeX 的 SHA-256 为 `b746049f90df2879b2bc7d8bcd6f616291b25551f5a664cf5aec61529321ca34`。

本轮先读当前稿件、关键代码和原始 JSON，再分别审方法、证据和文献，没有采用旧自动审稿分数。核查了核心公共库、配对比较、真实群体与规模实验、paradox、治理及执行限制代码；联网查阅最接近的原始论文和官方文档。全部结论仅针对该快照。

实际执行的轻量诊断及原始输出保存在：

- [check_evidence.py](/data/lvvyongji1/GWAS/review-stage/independent_review_20261015/check_evidence.py)
- [evidence_checks.json](/data/lvvyongji1/GWAS/review-stage/independent_review_20261015/evidence_checks.json)

该脚本重算了 S10 的 25 种子均值与同一 bootstrap 置信区间，核对了缓存形状、原始群体标签与池化定义，测试了 CPU 常量分数 AUROC，并在新子进程中对当前审计 hook 做了合成数据诊断。没有重跑完整 PCA/GWAS、GPU 规模网格或真实安全 GWAS 工具链；因此本轮是源代码和归档证据审计，加有限运行验证。

共享事实如下。下表各行均来自同一 S10 配对实验，而非拼接不同实验的均值。

| 攻击 | 25 种子平均 AUROC | 能支持的结论 |
|---|---:|---|
| 未使用已知成员贡献的 AF | 0.832436 | 当前基础 AF 实现 |
| 原始研究载荷的 PC 分数 | 0.878575 | 高于上述 AF |
| 原始 AF + PC | 0.935391 | 对上述 AF 的增量为 0.102955 |
| 利用相同已知成员信息的 conditioned AF | 0.940730 | 25/25 种子超过 PC，平均也超过现有组合 |

对 `combined - af_only` 重算的 95% bootstrap CI 为 `[0.0964710020, 0.1098678221]`，与归档逐项一致。**这个结果成立，但它没有测量 conditioned AF + PC 相对 conditioned AF 的增量。** 数据见 [s10_paired.json](/data/lvvyongji1/GWAS/results/raw/s10_paired.json)，实现见 [s10_paired.py](/data/lvvyongji1/GWAS/code/s10_paired.py:55)。

## 审稿意见一：侧重贡献与先例

**总体判断：倾向拒稿。** 基因组数据共享、联邦分析和隐私测量研究者会关心这个问题。项目的价值在于比较具体发布对象的风险，而当前稿件尚未证明这种比较带来足够新的、可推广的认识。

值得保留的是：作者承认了 PCA 攻击先例；对成员、参考集和独立非成员做了明确拆分；没有隐藏 conditioned AF 更强、真实混合人群中 AF 更强等不利结果。这些提高了结论的可信度。

**主要问题 1：攻击成立所需的明文发布对象，与核心引用不一致。**

稿件用 SF-GWAS 支撑“联邦 PCA 会共享研究拟合的 loading matrix”，见 [main.tex](/data/lvvyongji1/GWAS/paper/main.tex:200)，类似表述还出现在摘要、贡献与限制中。但 [SF-GWAS 原论文](https://davidfroelicher.com/publication/sfgwas/sfgwas.pdf) Methods 的 PDF 第 8 页写明：联合 PCs 保持加密用于后续计算，最终集体解密并共享的是关联统计量。其[官方输出说明](https://github.com/hhcho/sfgwas#output-format)也只列出关联结果和变异筛选输出。

[SF-PCA 的 Step 8](https://arxiv.org/html/2304.00129v1)允许按应用需求，经共同授权后解密 PCs 或投影结果；它不意味着必须向分析员发布明文载荷。可选解密场景可以研究，但必须明确配置、对象和接收者。

这不推翻“明文载荷可以泄漏”的实验现象，却直接影响论文对真实风险的定位。目前证据更适合支持条件命题：**当研究拟合的载荷及所需预处理信息被授权明文释放时，存在相应成员推断风险。** 不能据此声称已经证明所引用默认安全工作流有该攻击面。最低要求是复现一个具体工具版本和发布合同，说明攻击者实际能得到什么。相关公式还需要区分载荷 `k×M`、特征值 `k` 和个体投影 `n×k`，并写清均值、尺度、SNP 顺序和 QC mask 的来源。

**主要问题 2：PCA 泄漏机制本身已有直接先例，遗传学近邻文献也不完整。**

[Zari 等，PSD 2022](https://link.springer.com/chapter/10.1007/978-3-031-13945-1_19)已研究从可见 principal components 做成员推断，并指出小训练样本时风险更强。稿件 [main.tex](/data/lvvyongji1/GWAS/paper/main.tex:219)已正确承认这一点。因此，样本内拟合方向会区分成员、独立重拟合后效应下降以及增加 PCs 后风险增强，都不能单独承担主要新颖性。

更贴近本领域的 [Dervishi 等，Bioinformatics / ISMB 2023](https://academic.oup.com/bioinformatics/article/39/Supplement_1/i168/7210481)已经比较协作基因组 PCA 输出的成员风险与 AF 汇总风险。该文用公共参考模型和噪声化的个体投影，和本文研究拟合载荷不同；这一区别正应成为相关工作对照，而不能笼统声称中间对象从未被测量。

剩余的潜在贡献是：在明确条件下，定位 PCA 与 AF 风险排序随样本数、有效标记数、保留维度和群体结构变化的规律。目前大量网格证明了若干操作点，尚未形成足够稳定的解释或预测。机器学习主会不要求每篇文章都有定理，但这里需要比“某些格子 PC 更强、某些格子 AF 更强”更有迁移价值的发现。

**主要问题 3：三方拆分是必要控制，方法学原创性尚未建立。**

重叠参考与非成员造成约 0.132 AUROC 膨胀，能对应 [s1d_bias.json](/data/lvvyongji1/GWAS/results/raw/s1d_bias.json)。这值得报告。但稿件所谓 common shortcut 没有指出哪些已发表方法实际采用这一错误。

[Visscher 与 Hill，PLOS Genetics 2009](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000628)早已分析 proband、研究、参考样本的抽样关系、有限参考集、LD 和群体差异。它没有给出本文完全相同的 self-reference 偏差，因而不能说已经覆盖本文一切；但作者需要精确说明新增内容。当前更稳妥的定位是良好实验控制及其量化示例。

**成立贡献所需的最低补充：**纠正真实发布对象的引用；补齐上述近邻文献的逐对象比较；选定一个实际发布合同；用公平攻击比较证明一个超出已知 PCA 过拟合现象的新结论。仅扩大种子数不足以改变本审稿判断。

## 审稿意见二：侧重统计与实验设计

**总体判断：弱拒稿至拒稿。** 关键配对数字可靠，机制对照有解释力；但“诊断增加多少风险”和“风险边界在哪里”这两个核心问题，还没有在相同信息条件和可识别设计下完整回答。

真实优点包括：study/reference/nonmember 分离；S10 的方向校准与测试标签分开；独立 basis、另一站点 basis、载荷数量梯度和 msprime 检查。这些支持研究拟合子空间的样本内效应。本轮未发现 S10 把测试标签直接用于主配对方向选择。

**主要问题 1：公平强基线比较缺了一格。**

PC 校准和组合攻击知道目标研究约一半校准个体的成员标签。AF 使用同样已知成员基因型后，可以从发布频率中合法扣除其贡献，达到 0.940730。作者在 [main.tex](/data/lvvyongji1/GWAS/paper/main.tex:377)已经准确披露；这不是测试标签泄漏，也不应被批评成未承认错误。

真正缺失的是 `conditioned AF + 原始 PC` 对 `conditioned AF` 的增量。现有组合使用未条件化 AF，不能替代这一比较。真实数据也没有对应 conditioned AF 实验。约一半成员已知还是很强的辅助信息，应补少量已知成员比例，以及用独立 shadow study 校准的设置。各设置必须共享样本、SNP、预算和评价集。

**主要问题 2：引言询问最终关联结果之外的泄漏，核心比较实际上是 AF 分数。**

引言 [main.tex](/data/lvvyongji1/GWAS/paper/main.tex:129)询问中间诊断是否增加了最终统计量之外的信息；核心 final-only 在 [s2b_intermediate.py](/data/lvvyongji1/GWAS/code/s2b_intermediate.py:23)主要是 AF-LRT。真实数据无 phenotype，top-200 用 AF 差异代理，稿件 [main.tex](/data/lvvyongji1/GWAS/paper/main.tex:564)也已说明。

如果研究问题是 AF 发布与载荷发布的比较，应据此收窄标题和引言。如果要回答真实 GWAS 最终输出之外的新增风险，就需要针对实际输出的候选依赖攻击，再比较加入诊断前后。具体可用信息取决于真实工具，可能包括效应量、SE、p、N 或关联分数；不能把整个关联表对所有候选都相同，当作它无法与候选基因型联合产生攻击分数的证据。

**主要问题 3：sanitized-diagnostics paradox 的实现不支持目前的行为因果解读。**

[s25_paradox.py](/data/lvvyongji1/GWAS/code/s25_paradox.py:66)的 `decimals` 只映射到 PCA 保留维度，没有按相应小数位舍入谱；sanitized 分支直接返回停止。曝光却在 [同文件](/data/lvvyongji1/GWAS/code/s25_paradox.py:105)另行按同名小数位计算 AF-LRT。最终把 PCA 停止编号映射为 AF 精度，形成两条人为绑定曲线。质量平价只检查平均差值超过 -0.03，参数 `tol` 没有执行所暗示的等效性检验。

因此，当前实验可以描述一个构造的精度差示例，不能证明提供真实 sanitized 诊断会使分析者少请求危险输出。最低修改是降格表述；若保留主要贡献，就必须定义每一轮实际可见对象、统一继续/停止规则，并攻击真实累计发布内容。此问题不影响主 PC 泄漏现象。

**主要问题 4：“零结构”来自标签定义，不能当作遗传无结构实验。**

[s51_real_lowstruct.py](/data/lvvyongji1/GWAS/code/s51_real_lowstruct.py:65)把一个 super-population 中全部人标为 `pool`；[fst_like](/data/lvvyongji1/GWAS/code/s48_real_structure.py:136)在只有一个标签时必有组内、总体频率相同，返回零。这不排除真实亚群结构。本轮缓存核查显示 EUR/EAS/AFR 分别仍包含 5/5/7 个原始 population。

稿件 [main.tex](/data/lvvyongji1/GWAS/paper/main.tex:1030)承认 zero by construction，却随后沿用 largest real structure-zero cohorts。应改为“同一 super-population 的池化队列，组间统计量在当前标签定义下为零”，删除对真正零结构条件的外推。不能仅把参与者标签合并，就认为完成了零结构控制。

**主要问题 5：真实规模与结构边界仍是探索性结果。**

12 个单 population 队列确实有 PC>AF，但研究集仅 28–36 人，保留 10 PCs，接近几何饱和；这些事实也被作者披露。较大研究集的可信点包括 AFR n=350、reference=150、nonmember=161。EUR/EAS 的 n=350 各只有 6 个非成员，随后还分校准与评价，不能视为多个稳定复现。详见 [s56_real_poolmax.py](/data/lvvyongji1/GWAS/code/s56_real_poolmax.py:43)和配套 JSON。

固定每组人数没有固定总人数；另一组固定组成实验则让 study、reference 和评价集一起增长。五个拆分中 crossing 也不稳定：三次单一跨越、一个种子三次跨越、一个不跨越，见 [s60_share_crossing.json](/data/lvvyongji1/GWAS/results/raw/s60_share_crossing.json)。模拟 FST 边界还使用评估集上的 `max(AUC,1-AUC)` 选方向，见 [s21_sensitivity.py](/data/lvvyongji1/GWAS/code/s21_sensitivity.py:37)；作者已披露这种上偏，但摘要仍突出量化边界。

这些结果支持“排序依赖操作点”，不足以把 0.03–0.15 当作已确认且可推广的结构阈值。最低关键实验应固定 reference 大小、评价集大小和总研究人数，对结构与 n 分别做控制，并对 crossing 用独立方向校准。

**次要但应修正的问题：**同一 1000G 池的五次拆分 CI 描述的是该池上的拆分不确定性；同一参与者的不同染色体属于区域复现，不是独立人群研究。真实缓存部分实验只覆盖染色体片段，例如本轮核对的 chr22 为 16,051,249–23,512,554。GPU AUROC 没有正确处理 ties，已有单点无 ties 检查不能覆盖全网格；但也不能据此推翻所有已报告 GPU 数字。TPR@1% 是经验 ROC 点，若要给实际判定阈值，需要独立校准并报告实际误报区间。

## 审稿意见三：侧重工作流、治理和论文表达

**总体判断：倾向拒稿。** 发布清单、精度限制和可审计记录有实用意义；当前原型能演示部分约束逻辑，但尚不足以支撑完整治理机制和唯一出口的表述。

已有工作的优点是：manifest、egress monitor、hash chain 均有代码；作者明确不声称实现 MPC；有普通写文件、UDP、ctypes 的注入记录；也承认 bundle 的对象并非全部来自同一贡献者集合。这些是有价值的工程证据。

**主要问题 1：审计 hook 没有建立完整的唯一出口。**

[enforce.py](/data/lvvyongji1/GWAS/code/enforce.py:39)主要拦截 write-open、部分 socket/subprocess/ctypes 事件。[Python 官方文档](https://docs.python.org/3/library/sys.html#sys.addaudithook)明确指出，Python 层 audit hooks 不适合构造 sandbox；[PEP 578](https://peps.python.org/pep-0578/#why-not-a-sandbox)也区分可见性与隔离保证。

本轮提取当前 hook 的精确源码，在新子进程中预载 NumPy、安装 hook、丢弃 stdout/stderr，并用预开的 pipe 返回诊断。普通 `open(..., "wb")` 和 NumPy `tofile` 都被拦截；普通 `os.mkdir` 仍能在临时目录创建含合成值的目录名。该操作不需要 root、ctypes 或新原生代码，可形成不经过 monitor 的可观察文件系统变化。证据见 [evidence_checks.json](/data/lvvyongji1/GWAS/review-stage/independent_review_20261015/evidence_checks.json)。

这只是 exact-hook 层面的有限反例，**没有重跑整个研究运行时，也没有证明原始载荷数组已经从完整运行时外泄**。但足以否定“一条 pipe 是唯一可写渠道”的宽泛说法。如果威胁模型只允许一组封闭、已审查的运算，就应明确列出这组操作；如果允许普通 Python 研究代码，就需要操作系统权限层的隔离，或把贡献降格为拦截已列举误操作的原型。

**主要问题 2：manifest 约束的是标签与精度，未建立语义正确性或组合安全。**

[governance.py](/data/lvvyongji1/GWAS/code/governance.py:121)按 quantity ID 查策略并舍入数值，但 semantic_type 不构成内容验证；recipient 和 condition 主要检查传入值是否与声明相等，并未构成独立接收方鉴权和真实条件求值。跨运行关系来自 producer-asserted metadata，稿件已承认，不能当作系统从内容自动发现了 refinement。

Hash chain 可检测相对于可信记录头的修改；若整条链及链头由同一写入者掌控，单纯重算整链不会被链内校验发现。应定义谁持有可信锚点、谁能改 manifest 和 ledger，以及记录只能追加由谁保障。当前可以支持“实现声明与记录逻辑”，不能支持全面的输出安全保证。本文也没有正式组合隐私保证，不能把 AUROC 下降写成预算或安全证明。

**主要问题 3：风险降低没有与具体 GWAS 效用配对。**

删除载荷或舍入 AF 会限制某些攻击，工程上可以做到。但载荷也用于投影、校正与跨站点一致性。当前治理实验主要验证输出被阻断/舍入，以及现有攻击分数降低，没有在同一真实发布合同上测量校正、功效、校准和协作能力的损失。Paradox 中构造的质量分数不能替代该实验。

最低实用验证是比较当前工作流、拒绝发布研究载荷、改用独立公共 basis，以及适当的受控/噪声发布方案，同时报告实际分析效用和攻击风险。如果仅治理对外发布、内部加密计算仍可完成校正，也应明确内外边界，并通过实际工作流验证。

**主要问题 4：文章承载太多弱相关支线，掩盖了最值得讨论的问题。**

编译稿摘要按空白分词约 818 词，跨到第二页。它同时包含主结果、反例、构造控制、多个后续实验、脚本名和大量限制。读者难以在第一次阅读时识别唯一的主结论。章节名 The Diagnostics Release Dominates 也比“在特定操作点高于未条件化 AF”更宽。

建议主文只保留明确发布合同、同信息条件下的风险比较、一个可靠边界/机制结论。raw haplotype 自匹配、leave-self-out 亲缘统计、batch missingness 和描述性 epsilon 各有不同威胁模型；应大幅移入附录或删除无法帮助主论点的支线。尤其不能用大量限制句代替补齐决定性实验。

**成立系统贡献所需的最低补充：**明确可信计算边界和可执行操作；修正唯一出口表述；提供真实发布策略的风险—效用配对结果；记录版本、sample IDs、variant IDs、拆分及候选评分，使外部人员可复核关键结论。

## 综合裁决与修改优先级

三个侧重点的共识是：**现象可信度高于主贡献和实际部署结论的成熟度。** 当前最有力的内容是“研究拟合的载荷方向在若干真实与模拟条件下携带成员信号”，以及“不同发布对象可能提供互补攻击特征”。最薄弱的是把这些结果连到真实默认工作流、最强 AF 基线和通用治理结论。

影响接收判断的优先级如下。

| 优先级 | 必须解决的问题 | 验收标准 |
|---|---|---|
| P0 | 真实发布对象与文献冲突 | 具体工具版本、配置、对象、接收者和攻击者知识一致；纠正 SF-GWAS 引用 |
| P0 | 公平基线之外的增量未知 | 同一评价集上报告 conditioned AF、原始 PC、conditioned AF+PC；覆盖模拟及真实数据 |
| P0 | 主问题大于现有实验 | 加入实际最终输出攻击，或把全文明确收窄为 AF 与载荷发布比较 |
| P1 | 结构边界与零结构定义 | 去掉标签归零的外推；关键格子采用 held-out 方向；固定 n_ref、评价规模及总 n |
| P1 | paradox 与唯一出口过度解释 | 降格/移除现有结论，或补可见对象一致的反馈实验及实际隔离证据 |
| P1 | 缺少风险—效用权衡 | 对同一真实发布合同测量隐私攻击和 GWAS 分析效用 |
| P2 | 重现与表达 | 修复 GPU ties；保留关键评分和拆分；补近邻文献；重写摘要和贡献 |

可以先集中完成三个具有裁决力的实验，而不继续铺开支线：

1. **实际发布合同实验。** 从一个明确版本的工具出发，记录默认与可选发布配置。默认不提供明文载荷时，应承认现有 PC 攻击不适用；可选配置提供时，再测量它带来的实际风险和效用。
2. **同辅助信息的强基线实验。** 固定同一研究/参考/非成员和 SNP 集，对多个已知成员比例，比较 conditioned AF 与 conditioned AF+原始载荷。方向、模型和阈值用独立数据校准；保留逐候选分数。
3. **可识别的边界实验。** 固定 reference 与评价集，以固定总 n 改变结构，再固定组成改变 n。若主攻机器学习主会，进一步提出与 effective rank、k/n 或有效独立标记数有关的解释，并在未用于构建解释的数据上作预测验证。

与现有积累最接近的定位是具体工作流的基因组隐私测量。它仍需真实攻击面、实际危害与效用证据。若坚持机器学习主会，需要让经验比较导出可迁移认识；不必强行发明一个已有的 PCA 攻击。

文章可以收敛到这样的条件主张：**在明确授权释放研究拟合载荷的工作流中，载荷可能在 AF 之外增加成员信息；增量的大小与排序需在匹配攻击者知识、样本规模及群体组成后测量。** 其中“在强 AF 之外增加多少”仍是待完成项，应由新的公平实验决定，而不是先写为结论。

## 复核入口

从项目根目录运行以下轻量诊断；它只重写本独立审稿目录内的结果，临时合成文件会自动清理：

```bash
/opt/miniforge/bin/python -c 'import runpy; runpy.run_path("review-stage/independent_review_20261015/check_evidence.py", run_name="__main__")'
```

归档均值与 CI 的一致性验证已经通过。目录创建反例也已观测到。它们各自证明数值可核对和 hook 覆盖不完整，不能替代整套实验复现、真实工作流验证或普遍隐私保证。
