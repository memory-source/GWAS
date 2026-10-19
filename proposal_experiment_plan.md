# 面向 MPC 支撑的多站点 GWAS 工作流的可审计揭示治理：强制原型与经验泄露研究

> **Claim / Type**：PoPETs 系统原型 + 实证隐私研究。多站点 GWAS 缺的不是更快的安全计算，而是一层可预登记、可强制、可跨运行追踪的揭示治理——我们实现它，并测量它买到了什么、买不到什么。本工作实例化已被提出的去分类、审计与 provenance 思想，不提出新的密码学机制或新的隐私定义。
> **Problem**：MPC 只保证"除协议输出之外无泄露"，不决定输出集合本身是否安全；多站点 GWAS 的公开集合目前由实现细节和运行期临时判断隐式决定，无法事前承诺、无法强制、无法跨运行审计。
> **Fix**：运行前把公开集合固定为承诺清单（计算图 + 揭示策略 + 预承诺决策树），用受限 DSL 与中介运行时把它变成唯一可用出口，用只追加账本做跨运行审计，并对公开策略报告指定攻击模型下的经验泄露。
> **Contributions**：发布策略 × 揭示策略的效用–泄露经验前沿（核心实证） · 受限执行下的揭示强制原型与绕过注入评估 · 跨运行释放关系检测 · 与四类既有体系的逐能力对照。
> **Scope**：chr22 或 5k–20k 经 LD 剪枝的位点；使用公开人类基因组资源与模拟数据，不使用受控访问临床队列或新收集的受试者数据。

---

## 1. Idea 概述

### 1.1 问题动机

GWAS 的统计功效来自样本量，但基因型数据因法规（GDPR、我国人类遗传资源管理条例等）无法跨机构汇池。工程上成熟的做法是"把计算搬到数据所在处"：各站点本地计算汇总统计，只把汇总统计送出做 meta 分析；GA4GH TES 体系与欧盟 GDI 项目已在七个国家节点上跑通这条路。

这条路有隐含前提：假定汇总统计可以安全地公开。该前提在小站点下不成立——数百人规模的等位基因频率本身已足以反推个体。Pascoal 等人的系列工作量化了这一风险：某些配置下最多 8% 的基因组面临暴露风险；当多个研究数据集存在个体重叠时，最多 28.6% 的变异可被推断，其中 92.3% 会进一步使成员身份可判定。基因组信标服务同样是成熟的攻击面，针对它的成员推断与防御已有专门研究。

多方安全计算（MPC）降低了对参与方保守各自输入的信任要求——它并不删除全部信任，协议实现、运行时与执行环境仍在可信基之内。而 MPC 的模拟安全性只保证"除协议输出之外无泄露"——它不判断输出本身是否过度。现有隐私保护 GWAS 文献全部停在"算"这一层：从 Ostrak 等人把 EIGENSTRAT、FastPCA、EMMAX 与 genomic control 搬进 MPC 与 TEE，到 Aherrahrou 等人的分类综述，再到 GDI 的联邦流水线，处理的都是"这个统计量怎么在密文里算出来"。

已有工作分别研究了安全发布（Pascoal 系列）、语言级去分类、查询审计与工作流 provenance，但它们各自停在不同的层面。尚缺少一个系统把"运行前预登记揭示集合、运行时强制揭示出口、跨运行审计、基因组攻击评估"这四件事在 MPC 支撑的 GWAS 工作流上连接起来。本工作要回答的，正是这四件事连起来之后才出现的一组问题：一条 GWAS 工作流在运行过程中会公开哪些量？这些公开是否被事前固定？运行中实际发生的公开是否落在事前预登记的范围内？跨多次运行累积起来又会怎样？

这个问题在真实流程里更尖锐，因为 GWAS 里最有后果的动作不是"算"，是"判断"。群体分层该校正到几层主成分——加少了大量假阳性，加多了真实信号被当成分层抹掉。做这个判断所依据的量（λ_GC、主成分-表型相关、QQ 图形态）只有合池后才存在，任何单一站点都看不见。当前实践里，这类量的公开散落在实现代码中，由分析人员在运行期临时决定。

最后，用"预算"这个词本身就是危险的。Wagner 对基因组隐私度量的评估工作已经指出：多种被广泛使用的基因组隐私度量在对手强度增加时不满足单调性，因此不能作为可信的安全声明。这决定了本工作的基本立场——我们不声称提出一个隐私保证，我们提出一种可审计的承诺机制，外加一套在指定攻击模型下测得的经验性泄露估计。

### 1.2 核心想法

把研究单元从"算法"上移到公开集合的承诺与记账。

一条 MPC 化的 GWAS 工作流，在运行前被编译为一份承诺清单，它包含三部分：

1. 计算图：哪些步骤在明文、哪些在 MPC、每个节点的输出去向；
2. 揭示策略：每一项将被公开的量，连同精度、接收方、触发条件；
3. 预承诺决策树：所有分支事前写死，每条分支绑定其触发的公开子集。当分支条件依赖秘密值时，编译期只能得到所有可能路径的公开集合及其上界；某次运行实际走了哪条路径，需要运行后依据被记录的 condition bit 判定。

这份清单不是密码学意义上的 hiding commitment，而是一份带内容哈希、版本号与时间戳的已登记揭示清单（pre-registered reveal manifest）。它不隐藏任何东西，它做的事情是让"将要公开什么"在运行前就被固定下来，并且事后无法否认。

运行时，对由承诺清单生成的 DSL 工件而言，所有揭示必须经过一个统一的 reveal 原语；未经申报的揭示使工件通不过校验，无法进入执行。手工编写的任意代码不在这一保证范围内，1.3 的威胁模型会显式给出这条边界。运行后，审计层把实际公开序列与承诺比对，并把它写入一个跨运行的累计账本——因为真实的对抗场景不是一次干净的运行，而是同一个分析者在多轮之间用前一轮的输出去构造下一轮的策略。

这里必须区分两个层次，本工作的全部声明都挂在这个区分上：

- 句法层（可证明）：实际公开集合是否为承诺集合的子集。这是执行一致性，可以做到机械判定。
- 经验层（可测量，不可证明）：承诺集合在指定对手先验与攻击实现下的实际泄露量。这是估计值。

Agent 是可选的编译前端，不是这套设计存在的理由。承诺清单本身可以手工撰写、可以用模板生成、也可以由 LLM 从自然语言 SOP 编译得到。Agent 的价值主张被收缩为一句可检验的话：在把自然语言分析意图转成揭示策略的任务上，它能否比一个强规则/模板基线犯更少的安全关键错误。若实验显示不能，Agent 退为可选接口并被移出正文，不影响论文主贡献——这一点在第 3 节与第 9 节中都有对应的降级路线。

为什么以 MPC 作为执行底座？这里提出的是一个可检验的假设，而不是一条设计理由：MPC 的 reveal / open 是数量有限、位置明确的显式原语，因此"揭示"可以成为唯一需要被中介的出口；而普通联邦流程的外部输出散落在任意通信代码中，即便能通过 workflow engine 与 provenance hook 做记录，也不能把输出收敛到一个可控的出口集合。该假设的成立条件是清晰的：只有当 reveal / open 确实是唯一的外部输出、且 I/O 被中介运行时控制时，narrow enforcement point 这个说法才成立。我们把它写成假设，并且用 P0-5 的两个对照基线（logging-only provenance，以及一个用 Nextflow/CWL 风格 hook 插桩的联邦 pipeline）在同一组绕过注入下实测"阻断率、检测率、漏检率、TCB 规模"。若三条线在本实验中的表现接近，本工作就放弃"MPC 更适合"的说法，改为"MPC 是一个清晰的实例化后端"，把贡献收回到领域语义清单、攻击绑定评估与跨运行释放关系检测上。

### 1.3 方法概述

对手模型。我们显式区分三类对手，并说明系统对每一类防什么、不防什么。

- A1｜诚实但好奇的参与站点：遵循协议，试图从协议记录中推断他人输入。交给 MPC 的模拟安全性处理，本工作不额外贡献。
- A2｜受框架约束的工作流作者：漏写揭示、写错精度、把中间结果写进日志之外的文件、扩大接收方、放宽条件。本工作直接针对这一档，但有一个必须讲清楚的边界——A2 不能执行任意代码，不能直接调用底层 reveal / open，不能自由地做外部 I/O。
- A3｜自适应分析者：合法地提交多个承诺清单，用前几轮的输出构造下一轮分析，通过精度细化、条件放宽、重复公开同一统计量来累积放大泄露。本工作只提供检测与报告，不提供形式化的组合保证。

明确不防的：能执行任意代码、能改写执行框架或篡改日志的作者（这一档在 A2 之上，超出威胁模型）；绕过受控接口的带外信道与侧信道；机构被法律强制协力的情形；以及承诺集合本身可能就是一个过度公开的坏策略。最后这一条是本工作与"任何 allowlist 系统"的分界：系统保证的是"说到做到"，不是"说的就是对的"。

TCB 与强制机制。"执行层不可被绕过"不是一句可以随手写下的假设，而是一组必须实现、并需要被逐条验证的机制：

1. 运行路径上的可执行工件只能由承诺清单生成，不允许手工编写的任意代码进入执行路径；
2. 提供受限 DSL，其类型系统不暴露底层 reveal / open，也不提供文件写入与网络 API；
3. 所有外部 I/O 由中介运行时代理，沙箱阻断带外文件写入与网络访问，违约尝试被记录而不是被静默忽略；
4. 执行日志采用 hash chaining（只追加、可检测篡改），承诺清单按内容哈希绑定并带版本号与时间戳。

我们把可信计算基（TCB）界定为：schema 校验器、受限 DSL 的编译器或代码生成器、中介运行时与沙箱、日志链。只有这些组件被信任；工作流作者写下的任何内容都不在 TCB 内。C1 的作用域因此被严格绑定在由承诺清单生成的 DSL 工件上：能执行任意 Python、能直接调用底层 SPU reveal/open、能自由做 I/O 的作者不在 A2 之内，本文不声称能约束他们——覆盖这一档需要的是另一个层级的信息流与沙箱机制，那是另一项工作。P0-4 会针对这条边界做对抗测试，逐个报告哪些违约被阻止、哪些只能被检测、哪些落在 TCB 之外，并把 TCB 的规模与组成写成一张公开可检验的清单。

形式定义。以下三个定义在全文反复使用，此处一次性给出。

定义 1（公开对象与揭示策略）。一次公开由八元组刻画

  o = (quantity_id, semantic_type, precision, recipient, condition, run_id, release_time, derivation)

其中 semantic_type 说明被公开量的语义角色（汇总统计 / 校准诊断 / 质控统计），derivation 记录它由计算图中哪些节点推导而来。揭示策略 P 是四元组的有限集合 {(quantity_id, precision, recipient, condition)}。

精度有四种可形式化的形式：有效数字位数、分箱边界、阈值布尔（如"p < 5e-8 是否成立"）、top-k 截断。条件可能依赖秘密值；此时条件是否触发本身就是一次公开，因此分支判定必须计入公开记录，而不是只计分支内部的输出。

定义 2（执行一致性）。设 P 为承诺策略，O 为一次运行的只追加公开日志。定义两种一致性：

- 允许子集一致性：每个 o ∈ O 都能匹配某个 (q, π, r, c) ∈ P。它允许"承诺了但没公开"。
- 精确路径一致性：O 恰好等于预承诺决策树在该次输入下规定应当公开的集合。它要求"该公开的都公开了"。

两者的安全含义不同，实验中分别报告，不混用。

定义 3（账本关系）。设 L 为跨运行累计账本，在其上定义四种关系用于跨运行检测：细化（同一 quantity 出现精度严格更高的新记录）、条件放宽（同一 quantity 的 condition 从条件触发变为无条件）、重复（同一 quantity 在同一精度上被再次公开）、接收方扩张（同一 quantity 的接收方集合严格变大）。若一次新运行的公开集合在 L 上被既有记录支配，账本给出显式标记。

架构。

① 编译层（离线）：意图或手工/模板输入 → 计算图 + 揭示策略 + 预承诺决策树 → 通过 schema 与图校验。校验不通过则不产生可执行工件。这一层可以是一个 LLM Agent，也可以是一个模板编译器；两者在实验中作为不同实现被比较，且都不接触任何私有统计量。

② 执行层（在线，MPC）：在 SecretFlow/SPU（备选 MP-SPDZ）上执行计算图。所有揭示调用一个统一封装的 reveal 原语，该原语负责写入只追加日志。未在策略中申报的揭示在编译期即被拒绝，因此运行期不存在"想揭示但没申报"的合法路径。

③ 审计层（在线 + 离线）：逐条比对实际公开与承诺（定义 2）；把本轮公开追加进跨运行账本并按定义 3 检测关系；输出审计报告与经验性攻击面估计。

与现有工作的区别。相对隐私保护 GWAS 算法工作，我们不竞争 MPC kernel 的效率，而是把他们的算法当作执行层后端；相对联邦 GWAS 工程体系，我们处理"汇总统计不可安全公开"的情形，并把揭示变成一等公民；相对安全发布工作，我们把发布风险评估从"结果之后"前移到"流程编译期"，并扩展到工作流的所有中间公开与多次运行。

### 1.4 理论依据

依据一：MPC 的安全性不覆盖输出选择。协议保证的是"除规定输出外无泄露"，把"规定输出是什么"留给了应用层。在多站点基因组场景里，这个留白恰好落在最敏感的地方——汇总统计本身。

依据二：揭示在真实实现中是隐式的，因而不可审计。当前实践中，"公开了什么"散落在通信代码、日志、中间文件与人工判断里。把它变成一个显式对象是审计的前提——不是因为它更安全，而是因为它让"不安全"变得可被看见。

依据三：句法与经验必须分开，否则会做出站不住的声明。一个量公开后泄露多少，取决于对手的先验、参考面板、样本重叠、LD 结构与精度。Wagner 的工作表明常用的基因组隐私度量在这些因素变化时不满足单调性。因此本工作把"公开了什么"（可静态计算，可证明）与"公开导致什么"（只能在指定攻击模型下估计）分成两层，并明确不把前者称为隐私保证。

依据四：强制边界必须被证明存在，而不只是被声明。任何"未申报揭示会被拦下"的说法都依赖于一组具体机制（受限 DSL、中介运行时、日志链）。本工作把这条依赖关系写成 TCB 清单，并用违约注入去实测它的边界——这比声称一个覆盖所有情形的定理更有用，也更诚实。

依据五：预承诺的代价必须被测量，而不是被假设。把判断从运行期移到编译期，会在统计上付出代价——这是一个可测的量，且它的正负号事先并不确定。如果代价过大，结论就应该是"预承诺不适合这种分析"，而不是把结论包装成成功。

边界声明。本工作不提供差分隐私式的、与数据无关的隐私保证；不声称账本给出了组合预算；不声称"任何未申报公开都能检测"（该性质只在 TCB 完整时成立）；不声称"发布数量或公开精度单调对应泄露"。这些边界会在论文中显式列出。

### 1.5 预期贡献

C1｜清单生成工件上的 trace 合规性（在给定 TCB 下）。C1 是一条 trace 合规性质，不是端到端安全定理，作用域被严格限定在由承诺清单生成的 DSL 工件上。在该工件内、且 TCB 完整（揭示只经由统一原语、日志链 tamper-evident、工件生成路径上不存在绕过中介的 I/O）时，一次运行记录的公开集合满足允许子集一致性，判定代价与承诺策略规模线性相关。它不覆盖手工编写的任意代码，不类比差分隐私组合定理，也不声称给出任何数据无关的隐私界；对沙箱只能声称覆盖已测试过的 I/O 与逃逸路径，不声称阻断所有隐蔽信道。它的价值在于把这条依赖关系写成 TCB 清单，并用绕过注入逐条实测边界。

C2｜强制原型。受限 DSL 与代码生成器、静态检查、中介运行时与沙箱、hash 链日志，加上揭示策略与计算图的 schema；强制机制本身是 C2 的实质内容，而不是附带实现细节。开源仓库包含可在 3 站点模拟数据上端到端跑通的 demo。

C3｜发布策略 × 揭示策略的经验前沿。这是本文的核心实证贡献，也是对"只公开最终 p 值就安全"这一直觉的直接检验。我们把设计空间拆成两维并做小规模完整扫描：最终发布策略 F1–F5（仅显著命中 / top-k / 完整 p 值 / p 值 + 效应量 + 标准误 / 完整汇总统计），中间揭示策略 M0–M3（无诊断量 / 仅 λ_GC / 主成分相关诊断 / 缺失率统计）。对每个组合测量效用（功效、假阳性率、λ_GC 校准度）与经验泄露（成员推断 TPR@低 FPR、基因型恢复率），并报告中间揭示在固定最终发布条件下的边际泄露，且至少在保守档（F1/F2）与宽松档（F3/F4）各测一次。若边际泄露不显著，论文转向"发布策略前沿加可审计工作流清单"，中间揭示部分作为负结果保留。

C4｜跨运行释放关系检测。账本检测的是释放关系模式，而不是泄露量。按定义 3 识别细化发布、条件放宽、同量重复发布与接收方扩张，并检测"多次发布的并集覆盖了某次单发布所能给出的更高精度"这一模式。针对 A3，账本在这些模式出现时给出显式标记与累计视图。它只做检测与报告，不阻断；它与攻击风险的对应关系由 P1-1 的经验实验给出，而不是由定义给出。

论文正文承诺的贡献只有 C1–C4，其中 C3 承载全篇重心。编译质量基准与 Agent 前端只作为附录中的可选组件存在：基准集 pilot 规模 30–50 条，比较规则/模板编译器、LLM 零样本、LLM + schema 校验修复与人工撰写，指标是安全关键字段错误率（漏报揭示、多报揭示、精度错误、接收方错误、条件错误、路径公开量低估）；只有当它显示出显著优势时，才会在修订稿中被考虑升级为独立贡献。这是一条事前声明的降级路线，而不是事后解释。

### 1.6 相关工作定位

相对隐私保护 GWAS 算法工作（Ostrak et al. 2021；综述 Aherrahrou et al. 2024；较早期的方法学综述 Wu et al. 2022）：他们解决"统计量如何在密文里计算"。我们不与其在 MPC kernel 上竞争，而是把这类算法作为执行层可插拔后端。

相对安全发布工作（Pascoal 的 DyPS / I-GWAS / GenDPR）：他们评估发布之后的结果是否安全，给出动态发布的安全条件与成员推断风险度量。我们处理的是流程在运行过程中产生的全部公开，把它们前移到编译期作为承诺，并扩展到跨运行累积。两者互补：他们的安全条件可以作为我们揭示策略的可选约束。本条差异成立与否，取决于 P0-4 能否显示中间公开带来最终发布评估之外的新信息——我们把它当作待验证的前提，而不是既定事实。

相对语言级信息流与去分类（Sabelfeld & Myers 的综述；Sabelfeld & Sands 的 declassification 维度与原则；Li & Zdancewic 的 relaxed noninterference；Cruz & Tanter；Micinski 等人；Prasad 等人；Gancher 等人）：揭示策略在结构上是 declassification policy 的一个实例，一致性判定对应标准的非干涉类性质。差异在作用域：类型系统针对单个程序内的变量级流，我们针对跨机构工作流的揭示对象级承诺，载体是 MPC 执行日志而非程序语义，并且额外绑定基因组攻击评估。论文必须正面承认这一继承关系，否则会显得在重新发明 declassification allowlist。

相对 MPC 程序的静态分析（这一方向的代表工作仍以预印本为主，本工作不把核心论证建立在其上）：代码层检查"这一行的揭示是否安全"，流程层检查"这次运行整体承诺了什么"。二者互补而非竞争：一个成熟的代码级检查器完全可以作为本工作编译层的前置检查，把 unsafe reveal 挡在清单生成之前。

相对差分隐私与查询审计体系（PINQ、Airavat 的隐私账本与自适应查询审计传统，以及 Rogers 等人的 privacy odometer 一系工作）：概念上最接近。三处关键差异需要在论文中讲清楚。第一，我们的对手不是查询者，而是一组互不信任的数据持有方加一个自适应分析者。第二，泄露度量不是数据无关的 ε，而是在特定基因组攻击模型下的经验值。第三，PINQ/Airavat 通过 sensitivity 与噪声给出一条可组合的上界，我们的账本只提供检测与报告，没有数据无关上界——这是刻意的选择而非疏漏，因为 DP 在高维稀疏基因型上极难同时保住 GWAS 所需的效用，这正是该领域主流走向联邦与密码学方案的原因。

相对 workflow provenance 体系（W3C PROV 的实体—活动—代理模型；Nextflow / CWL 一类的科学工作流引擎所附带的数据血缘记录）：这类系统回答"这个结果由哪些输入与步骤产生"，其价值在可复现性与调试。它们对揭示的处理是观察式的——先发生，再记录。本工作要求的是先声明、后被允许发生，因此即使把 provenance 记录做到完备，也无法在运行前拒绝一个未申报的揭示。

相对透明度日志与供应链 attestation（Certificate Transparency 的 append-only Merkle 日志与 RFC 6962 的日志证明；Trillian 一类的可验证日志服务；in-toto、SLSA 一类的构建与制品溯源）：这类工作把"日志不可被悄悄篡改"做成了工程标准，本工作直接借用这一机制来实现日志链。差异在于信任链要回答的问题不同：它们证明"某个事件曾经被记录且未被改写"，本工作还需要证明"该事件在发生之前就被授权"。前者是完整性，后者是执行的出口约束。

相对受控统计分析平台（DataSHIELD 一类的"把分析送到数据处"框架）：这类系统的保护手段是限制可提交的统计操作集合，因此与"限制可公开的对象集合"在动机上部分重叠。差异有两点：它们的控制粒度是分析操作，我们的控制粒度是揭示对象及其精度、接收方与条件；它们不把运行期的中间诊断量纳入承诺范围，也不跨运行追踪释放关系。

下面的对照表把上述差异写成可逐格检验的形式，而不是一句"我们更强"。

| 体系 | 记录了什么 | 强制唯一揭示出口 | 带工作流对象语义 | 绑定攻击模型评估 |
|---|---|---|---|---|
| workflow provenance / lineage | 数据与步骤的血缘 | 否 | 否 | 否 |
| 透明度日志 / 供应链 attestation | 事件与构建的可验证存在性 | 否 | 否 | 否 |
| 语言级 IFC / 去分类（Jif、FlowCaml、Jeeves 一系） | 程序内的变量级流 | 是（程序内部） | 否（无 GWAS 语义字段） | 否 |
| DP 查询账本与隐私里程计（PINQ、Airavat） | 查询与预算消耗 | 是（查询接口内） | 否 | 否 |
| MPC 代码级 unsafe-reveal 检查 | 代码行的揭示安全性 | 否（只做检查） | 否 | 否 |
| 本工作 | 跨机构工作流的揭示对象与释放关系 | 是（受限工件内） | 是 | 是 |

这张表同时说明了本工作不主张什么：我们不提出新的强制机制原理（来自 IFC 与查询审计），不提出新的日志完整性机制（来自透明度日志），也不提出新的攻击（来自基因组隐私文献）。我们主张的是把这四类已有能力在 MPC 支撑的 GWAS 工作流上连成一个可运行、可审计、可被经验检验的整体，并且用实测指出每一格的边界在哪里。P0-5 的 logging-only 基线与联邦插桩基线把前两行的"否"直接测出来。

相对 LLM × 密码学工作（TFHE-Coder；Privformer 给出的 MPC 下 Transformer 推理代价量级）：我们不把模型放进 MPC 内部。Privformer 报告的代价量级（LAN 环境下翻译 64 token 约 19 分钟）是这一设计决策的直接依据。

### 1.7 参考文献

> 以下条目均经检索工具核验标题、作者与出处存在。带 `[v]` 者在本次写作中已核验。

- [1] [v] A. Ostrak, J. Randmets, V. Sokk, S. Laur, L. Kamm, "Implementing Privacy-Preserving Genotype Analysis with Consideration for Population Stratification", *Cryptography*, 5(3):21, 2021. https://doi.org/10.3390/cryptography5030021
- [2] [v] N. Aherrahrou, H. Tairi, Z. Aherrahrou, "Genomic privacy preservation in genome-wide association studies: taxonomy, limitations, challenges, and vision", *Briefings in Bioinformatics*, 2024. https://doi.org/10.1093/bib/bbae356
- [3] [v] T. A. Pascoal, J. Decouchant, M. Völp, "Secure and distributed assessment of privacy-preserving GWAS releases (GenDPR)", *Middleware*, 2022. https://doi.org/10.1145/3528535.3565253
- [4] [v] T. Pascoal, "Secure, privacy-preserving and practical collaborative Genome-Wide Association Studies" (PhD thesis: DyPS / I-GWAS / GenDPR), University of Luxembourg, 2022. https://orbilu.uni.lu/handle/10993/51996
- [5] [v] I. Wagner, "Evaluating the Strength of Genomic Privacy Metrics", *ACM Transactions on Privacy and Security*, 2016. https://doi.org/10.1145/3020003
- [6] [v] R. Venkatesaramani, Z. Wan, B. A. Malin, Y. Vorobeychik, "Defending Against Membership Inference Attacks on Beacon Services", *ACM Transactions on Privacy and Security*, 2023. https://doi.org/10.1145/3603627
- [7] [v] Z. Fang, "LSTM for membership inference attacks on genomic AAF pools and beacon services", 2024. https://doi.org/10.17760/d20705206
- [8] [v] N. Carlini, S. Chien, M. Nasr, S. Song, A. Terzis, F. Tramèr, "Membership Inference Attacks From First Principles", *IEEE S&P*, 2022. https://doi.org/10.1109/sp46214.2022.9833649
- [9] [v] R. Cruz, É. Tanter, "Existential Types for Relaxed Noninterference", *APLAS*, 2019. https://doi.org/10.1007/978-3-030-34175-6_5
- [10] [v] K. K. Micinski, J. Fetter-Degges, J. Jeon, J. S. Foster, M. R. Clarkson, "Checking Interaction-Based Declassification Policies for Android Using Symbolic Execution", *ESORICS*, 2015. https://doi.org/10.1007/978-3-319-24177-7_26
- [11] [v] S. Prasad, R. Yerraguntla, S. Sharma, "Secure Information Flow Typing in LUSTRE", 2022. https://arxiv.org/abs/2201.00184
- [12] [v] J. Gancher, S. Gibson, P. Singh, S. Dharanikota, B. Parno, "Owl: Compositional Verification of Security Protocols via an Information-Flow Type System", *IEEE S&P*, 2023. https://doi.org/10.1109/sp46215.2023.10179477
- [13] [v] Y. Zhang, M. Zheng, Q. Lou, "MPC-Patch-Bench: Security-Aware LLM Code Patch for Multi-Party Computation", arXiv:2606.11416, 2026. https://arxiv.org/abs/2606.11416
- [14] [v] M. Kumar, J. Xue, M. Zheng, Q. Lou, "TFHE-Coder: Evaluating LLM-agentic Fully Homomorphic Encryption Code Generation", arXiv:2503.12217, 2025. https://arxiv.org/abs/2503.12217
- [15] [v] Y. Akimoto, K. Fukuchi, Y. Akimoto, J. Sakuma, "Privformer: Privacy-preserving Transformer with MPC", *IEEE EuroS&P*, 2023. https://doi.org/10.1109/eurosp57164.2023.00031
- [16] [v] A. Chunduru, M. Joel, Z. Li, R. Madduri, "Towards Global Federated Genome-Wide Association Meta-Analysis Using GA4GH TES", arXiv:2609.02227, 2026. https://arxiv.org/abs/2609.02227
- [17] [v] C. Hernandez-Ferrer et al., "GDI D8.9 – Distributed analysis and federated learning PoC", Zenodo, 2025. https://doi.org/10.5281/zenodo.16882880
- [18] [v] C. Hernandez-Ferrer, N. Diaz, S. Capella, L. Portell, "GDI D8.1 – Federated access and federated learning technologies for advanced data analytics scenarios", Zenodo, 2025. https://doi.org/10.5281/zenodo.17989131
- [19] [v] B. C. Akdeniz et al., "COSGAP: COntainerized Statistical Genetics Analysis Pipelines", *Bioinformatics Advances*, 2024. https://doi.org/10.1093/bioadv/vbae067
- [20] [v] Z. Wu, W. Xiao, H. Deng, L. Yang, "Research on data privacy protection methods based on genome-wide association study", *Int. J. Embedded Systems*, 2022. https://doi.org/10.1504/ijes.2022.122105
- [21] [v] S. Zapechnikov, "Secure multi-party computations for privacy-preserving machine learning", `BICA*AI`, 2022. https://doi.org/10.1016/j.procs.2022.11.100
- [22] [v] R. M. Rogers, S. P. Vadhan, A. Roth, J. Ullman, "Privacy Odometers and Filters: Pay-as-you-Go Composition", *NeurIPS*, 2016. https://arxiv.org/abs/1605.08294
- [23] [v] N. Homer, S. Szelinger, M. Redman, et al., "Resolving individuals contributing trace amounts of DNA to highly complex mixtures using high-density SNP genotyping microarrays", *PLoS Genetics*, 4(8):e1000167, 2008. https://doi.org/10.1371/journal.pgen.1000167
- [24] [v] S. Sankararaman, G. Obozinski, M. I. Jordan, E. Halperin, "Genomic privacy and limits of individual detection in a pool", *Nature Genetics*, 41(9):965–967, 2009. https://doi.org/10.1038/ng.436
- [25] [v] S. S. Shringarpure, C. D. Bustamante, "Privacy Risks from Genomic Data-Sharing Beacons", *American Journal of Human Genetics*, 97(5):631–646, 2015. https://doi.org/10.1016/j.ajhg.2015.09.010
- [26] [v] J. L. Raisaro, F. Tramèr, Z. Ji, et al., "Addressing Beacon re-identification attacks: quantification and mitigation of privacy risks", *JAMIA*, 24(4):799–805, 2017. https://doi.org/10.1093/jamia/ocw167
- [27] [v] A. Sabelfeld, A. C. Myers, "Language-Based Information-Flow Security", *IEEE Journal on Selected Areas in Communications*, 21(1):5–19, 2003. https://doi.org/10.1109/JSAC.2002.806121
- [28] [v] A. Sabelfeld, D. Sands, "Dimensions and Principles of Declassification", *IEEE CSFW*, 2005. https://doi.org/10.1109/CSFW.2005.15
- [29] [v] P. Li, S. Zdancewic, "Downgrading policies and relaxed noninterference", *POPL*, 2005. https://doi.org/10.1145/1040305.1040319
- [30] [v] F. D. McSherry, "Privacy integrated queries: an extensible platform for privacy-preserving data analysis", *SIGMOD*, 2009. https://doi.org/10.1145/1559845.1559850
- [31] [v] I. Roy, S. T. V. Setty, A. Kilzer, V. Shmatikov, E. Witchel, "Airavat: Security and Privacy for MapReduce", *NSDI*, 2010. https://www.usenix.org/legacy/events/nsdi10/tech/full_papers/roy.pdf

- [32] W3C, "PROV-O: The PROV Ontology", W3C Recommendation, 2013. https://www.w3.org/TR/prov-o/
- [33] P. Di Tommaso, M. Chatzou, E. W. Floden, P. P. Barja, E. Palumbo, C. Notredame, "Nextflow enables reproducible computational workflows", *Nature Biotechnology*, 35:316–319, 2017. https://doi.org/10.1038/nbt.3820
- [34] B. Laurie, A. Langley, E. Kasper, "Certificate Transparency", RFC 6962, IETF, 2013. https://www.rfc-editor.org/rfc/rfc6962
- [35] J. Torres-Arias, H. Afzali, T. K. Kuppusamy, R. Curtmola, J. Cappos, "in-toto: Providing farm-to-table guarantees for bits and bytes", *USENIX Security*, 2019. https://www.usenix.org/conference/usenixsecurity19/presentation/torres-arias
- [36] B. Gaye, Y. Marcon, J. Isaeva, et al., "DataSHIELD: taking the analysis to the data, not the data to the analysis", *International Journal of Epidemiology*, 43(6):1929–1944, 2014. https://doi.org/10.1093/ije/dyu188
- [37] OpenDP, "The OpenDP Project", 2020–. https://opendp.org/

**说明**：文献 [13]、[14]、[16]、[17]、[18] 属于近一两年出现、以预印本或项目交付物形式存在的工作，本工作只把其中与本主题直接相关的部分作为领域现状的旁证，核心论证不建立在其上。文献 [32]–[37] 是标准参考，正式写作时会逐条核对版次与页码。Pascoal 系列中安全发布条件的形式定义细节仍需在撰写相关工作与攻击设定时逐条对齐。

---

## 2. 复杂度评估

- 综合难度：⭐⭐⭐⭐ — 难点不在算法，而在四件事必须同时成立：MPC 执行、强制机制、攻击实现、统计基线。任何一块不可信都会污染结论。
- 算力：不依赖大规模 GPU。建议 ≥ 16 核 CPU、≥ 64 GB 内存、≥ 500 GB 存储；4 站点可在单机多进程模拟。
- 数据：低难度。1000 Genomes 与 msprime 均公开。不使用受控访问的临床队列或新收集的受试者数据——这是伦理与合规上的实质优势，应在论文首页声明。
- 实现：大。约 55% 可复用（PLINK2、SPU/MP-SPDZ、统计库），45% 需从零写（受限 DSL 与校验器、中介运行时、审计器与账本、攻击实现、基准集）。
- 前置知识：需补 GWAS 流程（QC → PCA → 关联检验 → 多重检验校正）与群体分层、λ_GC 的含义；需掌握一套 MPC 框架（建议先 SPU）。已有隐私计算基础是当前最大优势。

scope 纪律：以 chr22 或 5k–20k 经 LD 剪枝的位点为规模上限，定位为隐私与工作流案例研究，不追求生产级 GWAS 引擎。

---

## 3. 实验目标（按优先级排序）

优先级分三段。P0-1 → P0-2 → P0-3 → P0-4 决定项目是否值得继续，最先跑；P0-5 与 P0-6 与它们并行，用来证明系统本身成立。P0-7 与 P0-8 是两项主贡献实验的最小版本，只在 P0-4 有信号后才需要完整展开。

### P0（核心，必须完成）

- **P0-1 — 明文金标准与站点切分**
  - Goal: 建立后续一切对比的效用参照物，并构造隐私风险主要来源的小站点。
  - Setup: 1000 Genomes 三个祖先面板切为 5 个站点（两个纯 EAS、两个纯 EUR、一个 AFR 或混血），每站点 200–800 人；植入一个主效位点与若干小效应位点，分层强度作为可控旋钮；chr22 或 LD 剪枝后的 5k–20k 位点。用 PLINK2 在合池数据上跑完整流程。
  - Pass: 合池 λ_GC 落在 [0.95, 1.10]，植入位点可被检出，功效曲线随效应量单调。
  - Fail: 检查表型模拟与分层参数；若金标准本身校准不良，后续所有效用结论作废。

- **P0-2 — 联邦 meta 四档基线的效用与泄露**
  - Goal: 建立"各算各的再 meta"这条现实工程主流的效用与泄露基线，同时为 MPC 底座的必要性提供实测证据。
  - Setup: 不以 naive 逆方差加权作为唯一样本，至少四档：固定效应逆方差加权、随机效应 meta、按祖先分层 meta、带站点协变量或 genomic control 的 meta。在祖先构成不平衡的切分上运行，并记录每个站点发出的全部汇总统计作为公开集合。
  - Pass: 报告出在大分层不平衡条件下假阳性相对金标准的膨胀倍数（这是全篇动机的关键数字），且该膨胀在四档修正后仍有至少一档未被消除。
  - Fail: 若所有现实修正都能消除膨胀，说明普通联邦路线的效用问题不成立。此时 MPC 不再作为必要前提，动机改由泄露侧单独支撑，MPC 收窄为"一个特别适合实现揭示治理的后端"。

- **P0-3 — 攻击端点校验与祖先匹配校准**
  - Goal: 建立可信的泄露度量，在先打穿"完整汇总统计"这一端点的同时，排除"攻击其实是祖先判别"这一混淆。
  - Setup: 至少两类攻击族并报告：(i) Homer / Sankararaman 风格的聚合等位基因频率成员推断；(ii) 信标/AAF 路线上的似然比或机器学习攻击。报告 TPR@FPR≤0.1% 与 ≤1% 两档，不使用平均准确率。对手先验覆盖强弱两档。校准子项：member/non-member 必须祖先匹配，参考面板重叠与不重叠各报一档，逐祖先分别报告 TPR，位点按 LD 剪枝或 LD block 处理，并说明攻击者是否知道表型/病例对照状态。基因型恢复全文只用一个主定义。
  - Pass: 在 F5 端点上至少一类攻击在低 FPR 下显著优于随机；强先验明显强于弱先验；在祖先匹配、面板不重叠的严格控制下攻击增益仍然存在；并且一个只使用祖先标签或 PCA 投影的对照攻击被证明明显弱于真实攻击——后者作为混淆上界，用来排除"攻击增益其实就是祖先判别"。
  - Fail: 若攻击在全公开上无效，或增益在祖先匹配后消失，则不允许进入任何隐私结论。先修攻击、对齐 Homer/Sankararaman 或 GenDPR 的设定，再重做。

- **P0-4 — 中间揭示的边际泄露**
  - Goal: 直接检验"只公开最终 p 值与效应量即为低泄露"这一直觉。它决定论文的主叙事方向，必须最早给出答案。
  - Setup: 固定最终发布策略，逐档加入中间揭示——仅 λ_GC、主成分特征值谱与主成分-表型相关、各位点缺失率与 QC 统计、以及 secret-dependent 分支的 condition bit。每档运行 P0-3 的攻击并报告相对前一档的增量。至少在两类最终发布档下各测一遍：保守档（F1 或 F2）与宽松档（F3 或 F4）。
  - Pass: 在至少一档最终发布下，中间揭示的增加带来方向一致、可测的边际泄露提升，且至少一类语义不同的诊断量有显著边际效应。
  - Fail: 若最终发布主导泄露而中间揭示的边际效应可忽略，论文立即重构为"多站点 GWAS 工作流的发布策略经验前沿"——比较 F1–F5 的隐私风险与效用；中间揭示作为负结果保留（常见诊断量在当前攻击模型下边际风险有限，但系统仍保证它们被记录与审计），揭示策略机制保留但在标题与贡献中的权重下调。

- **P0-5 — 强制边界与绕过注入**
  - Goal: 验证系统把揭示变成受限工件内的唯一出口，并明确哪些行为被阻止、哪些只被检测、哪些落在 TCB 之外。
  - Setup: 分两层注入。策略层：未申报揭示、提高精度、更换接收方、条件不匹配、跨运行重复揭示。实现层（在真实运行时上执行，而非仅 mock）：直接调用底层 SPU reveal/open、写临时文件、写 stdout/stderr 或异常信息、发起网络请求、启用 subprocess、利用序列化副作用、尝试导入逃逸、篡改本地日志后再提交。同时实现两个对照：logging-only provenance 基线，以及一个用 Nextflow/CWL 风格 hook 插桩的普通联邦 pipeline。
  - Pass: TCB 覆盖范围内的违约 100% 被阻止或检测，被阻止与被检测的比例分别报告；logging-only 与联邦插桩基线在"拦截"维度上为零；TCB 之外的行为被显式列出而非回避；同时报告 TCB 规模（组件清单与代码行数）。
  - Fail: 先修静态检查、运行时中介与日志链。若无法防住任意代码作者，就把威胁模型正式收窄为"只能提交由清单生成的 DSL 工件"，把这一条写进威胁模型而不是含糊带过。

- **P0-6 — MPC 揭示绑定演示**
  - Goal: 证明承诺清单与 reveal 原语能绑定到真实的 MPC 后端，而不是只在 mock 运行时上成立。
  - Setup: 2–3 站点、每站点约 200 人、1k–5k 位点，在 SecretFlow/SPU 上只实现一个可核验的最小目标——安全聚合等位基因计数并据此算出逐位点的关联统计量与 p 值。不做安全主成分，不做完整 EMMAX。揭示严格经由统一原语。
  - Pass: 与明文金标准数值高度一致（逐位点 p 值相关系数 > 0.99），所有揭示进入审计日志，策略在真实执行中生效。
  - Fail: 若该 demo 超预算，主线退回明文模拟的强制边界演示，并在论文中明确说明 MPC 后端只是绑定证明点而非性能贡献。

- **P0-7 — 跨运行释放关系的最小检测**
  - Goal: 证明账本不是 append-only provenance，而能判定跨运行的释放关系。这是 C4 能否作为主贡献的关键证据，只在最小规模上做。
  - Setup: 3 轮自适应承诺清单，分别覆盖精度细化、条件放宽、接收方扩张三类关系。按定义 3 判定，并至少在其中一类关系上测量攻击风险随轮次的变化。
  - Pass: 三类关系均被正确标记，且至少一类关系的标记方向与攻击风险变化方向一致；关掉账本后这些模式不可见。
  - Fail: 若账本只能记录不能判定，或标记与风险完全无关，C4 降级为 provenance 特性并从主贡献中移除。

- **P0-8 — 策略字段的最小消融**
  - Goal: 证明 precision 与 condition 不是 schema 装饰，而是必须被记录与审计的一等对象。
  - Setup: 只选一个高风险公开量（等位基因频率或 p 值加效应量），比较 full precision、四舍五入、分箱、阈值布尔、top-k 五种形式；另对同一公开量比较 conditional 与 always，以及 condition bit 可见与不可见。
  - Pass: 至少一种形式变化带来可测的攻击差异；若呈非单调但有可解释结构，同样算通过，只要该结构被重复观测到。
  - Fail: 若字段变化对攻击风险没有稳定影响，precision 一律降级为审计元数据，不得写入"精度是隐私调节旋钮"这类表述。

### P1（重要，建议完成）

- **P1-1 — 跨运行释放关系检测（扩展）**
  - Goal: 在 P0-7 的最小版之上，把账本判定扩展到完整关系集与更长轮次。
  - Setup: 5 轮以上的自适应承诺清单，覆盖四类关系加并集覆盖模式，并把账本判定与攻击风险曲线逐轮对齐。
  - Pass: 四类关系与并集覆盖均被正确标记，且标记与风险变化的对应关系在多数轮次上成立。
  - Fail: 回退到 P0-7 的最小版，只报告已验证的那一类关系。

- **P1-2 — 策略字段的完整消融**
  - Goal: 在 P0-8 的最小版之上，把精度与条件的语义检验扩展到多类公开量。
  - Setup: 对 p 值、效应量、标准误、等位基因频率、缺失率与诊断量分别做五档精度形式与两种条件模式的交叉，报告每类攻击下的变化。
  - Pass: 至少两类公开量呈现稳定可解释的风险—精度关系；若普遍非单调，则给出非单调性的结构解释。
  - Fail: 全部降级为审计元数据，论文明确写出"我们不把精度当作隐私控制手段"。

- **P1-3 — 预承诺的统计代价**
  - Goal: 回答"运行期不做临机判断，会不会毁掉 GWAS 效用"。
  - Setup: 三种流程对比——自适应分析者在看到诊断量后自行选择主成分数、预承诺决策树（主成分数取自事前声明的 {2,5,10} 及规则）、固定保守流程（恒定 10 个主成分）。至少两档祖先不平衡与两档效应量。指标为功效、假阳性率、λ_GC、显著发现数。
  - Pass: 预承诺决策树接近自适应流程，且明显优于固定保守流程。
  - Fail: 若代价过大，论文转向量化预承诺的代价与适用边界，讨论何种分析适合预承诺、何种需要噪声化诊断量或受限交互。负结果在此是贡献。

- **P1-4 — 发布策略 × 揭示策略二维扫描**
  - Goal: 得到效用–泄露–开销的经验前沿，验证揭示策略的粒度确实影响泄露。
  - Setup: 先做最小网格（2 档最终发布 × 4 档中间揭示）；若 P0-4 有信号，再扩到全交叉 F1–F5 × M0–M3。
  - Pass: 前沿单调或呈现可解释的非单调结构，且中间揭示在部分最终发布档下才可见边际风险。
  - Fail: 若所有合理发布策略都严重泄露，转向"安全发布策略设计"，而非"工作流中间揭示审计"。

- **P1-5 — CompileBench pilot**
  - Goal: 只在有数据的情况下回答"Agent 前端是否值得保留"。
  - Setup: 30–50 条高质量标注意图，比较规则/模板编译器、LLM 零样本、LLM + schema 校验修复、人工撰写。指标为安全关键字段错误率。
  - Pass: LLM + schema 在安全关键错误率上显著优于模板基线。
  - Fail: Agent 停留在附录或降为可选接口，主贡献保持 C1–C4 不变。

### P2（补充，时间允许时完成）

- P2-1 — 站点数与样本量敏感性：只取少数关键点（站点数 2/4、每站点 100/500 人、LD 剪枝强度 2 档），并至少在 1000G 真实祖先与 msprime 可控结构两套数据上确认趋势方向一致。
- P2-2 — CompileBench 全量扩展：仅当 P1-5 显示显著优势时才扩到 100–200 条。
- P2-3 — TEE 定性对照表；不实现。

### P3（可选）

- P3-1 — 审计报告一页可视化（本轮公开了什么、触发哪些账本关系）。
- P3-2 — 易被误判为无害的公开量清单，整理经 P0-4 与 P1-2 验证的失败模式。

---

## 4. 资源情况

### 可用代码与框架

- PLINK 2.0（https://www.cog-genomics.org/plink/2.0/）— 明文基线：QC、主成分、关联检验。公开，单文件可执行。
- SecretFlow / SPU（https://github.com/secretflow/spu）— 执行层首选，Python 友好、中文文档完善，学习曲线最平缓。Apache-2.0。
- MP-SPDZ（https://github.com/data61/MP-SPDZ）— 执行层备选，协议最全、性能最好，学习曲线陡。仅在 SPU 性能不足时启用。
- msprime（https://github.com/tskit-dev/msprime）— 模拟祖先结构与基因型。
- scikit-allel / pandas-plink / statsmodels — 统计与 IO。
- 参考实现：Ostrak et al. 2021（文献 [1]）的安全主成分与分层校正思路。检索未见其公开仓库，需自行实现，这是本项目最大的工程风险。

### 数据集

| 数据集 | 来源 | 规模 | 用途 | 获取难度 |
|---|---|---|---|---|
| 1000 Genomes Phase 3 | https://www.internationalgenome.org/data | 2504 个体、约 260 万 SNP | 祖先结构骨架与站点切分 | 公开 |
| msprime 模拟 | https://github.com/tskit-dev/msprime | 任意 | 可控祖先分层与信号植入 | 公开 |
| GWAS Catalog | https://www.ebi.ac.uk/gwas/ | 大量已发表汇总统计 | 挑选有已知信号的性状作为阳性对照 | 公开 |
| Neale Lab UKB 汇总统计 | http://www.nealelab.is/uk-biobank | 大规模 | 可选：更强信号来源 | 公开（仅汇总统计） |

### 预训练模型

| 模型 | 来源 | 用途 | 获取难度 |
|---|---|---|---|
| 商用 LLM API | API | 编译层与审计层的实现之一 | 公开（需 key） |

不做本地模型消融：那属于 LLM 代码生成方向，与本工作的 PETs 主线无关。

---

## 5. 环境配置

- 硬件：主实验机 ≥ 16 核 CPU、≥ 64 GB 内存、≥ 500 GB SSD。多机 MPC（4 台同配置、千兆局域网）推荐但非必需。不需要 A100/H100 集群。
- 软件：Python 3.10；统计侧 pandas / numpy / scipy / statsmodels / scikit-allel / pandas-plink / msprime / tskit；隐私侧 secretflow（含 SPU）；Agent 侧 openai / anthropic / pydantic（schema 校验）；沙箱与日志侧 Docker 或 bubblewrap 加 hash 链实现；评估侧 scikit-learn / matplotlib / seaborn。
- 安装步骤

```bash
conda create -n reveal-gov python=3.10 -y && conda activate reveal-gov
pip install pandas numpy scipy statsmodels scikit-allel pandas-plink msprime tskit
pip install -U secretflow
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip
unzip plink2_linux_x86_64_latest.zip && chmod +x plink2 && sudo mv plink2 /usr/local/bin/
pip install openai anthropic pydantic scikit-learn matplotlib seaborn jupyter
# 备选后端，仅在 SPU 性能不足时启用
# git clone https://github.com/data61/MP-SPDZ.git && cd MP-SPDZ && make -j
```

- 冒烟测试

```bash
plink2 --version
python -c "import msprime, secretflow; print(msprime.simulate(10, random_seed=1).num_sites, secretflow.__version__)"
plink2 --vcf chr22_subset.vcf.gz --maf 0.05 --geno 0.05 --make-bed --out siteA_qc
plink2 --bfile siteA_qc --pca 10 --out siteA_pca
```

---

## 6. 实验步骤

阶段划分与 §3 的 P0/P1 编号一一对应，此处只列工件与验收顺序，不重复实验设计。

### 阶段 1：环境与数据（3–4 周）
- [ ] 环境与冒烟测试；下载 1000 Genomes chr22 子集
- [ ] 补齐 GWAS 领域知识：读一篇标准流程综述并用 PLINK2 亲手跑通完整流程（硬性验收）
- [ ] `make_sites.py`、`plant_signal.py`；合池 λ_GC 与植入信号通过 P0-1 的 Pass 条件

### 阶段 2：核心实现（8–12 周）
- [ ] 2A（最先做，3 周）：`pipeline_plaintext.py`、四档 `federated_meta.py`、`attacks/membership.py` 与 `attacks/recovery.py`（含祖先混淆上界）→ 完成 P0-1 至 P0-3
- [ ] 2B（3 周）：`schema.py`、受限 DSL 与代码生成器、`validator.py`、`runtime/mediation.py` 与 `sandbox.py`、hash 链 `ledger.py`、logging-only 与联邦插桩两条对照基线 → 完成 P0-5
- [ ] 2C（2 周）：`auditor.py`（两种一致性分别判定）、`ledger_relations.py`（定义 3 四类关系与并集覆盖）、自适应分析者模拟器 → 完成 P0-7
- [ ] 2D（4–6 周，关键路径）：`pipeline_mpc.py`（SPU，安全聚合等位基因计数加逐位点关联统计量）、与金标准对齐 → 完成 P0-6；超预算即触发降级
- [ ] 2E（2 周，可并行）：`compiler.py`、`attack_surface.py`；3 站点 × 200 人 × 5k 位点端到端跑通

### 阶段 3：实验运行（8–10 周）
- [ ] P0-4（最先跑，决定叙事）；P0-8 最小消融
- [ ] P0-2、P0-3、P0-5、P0-6、P0-7 复跑，3 次独立重复
- [ ] P1-1 至 P1-4；P1-5 pilot

### 阶段 4：结果分析（4 周）
- [ ] 三张核心图：效用–揭示策略、泄露–揭示策略、开销–规模
- [ ] TCB 三列表：被阻止 / 被检测 / 不可防
- [ ] 判定 C3 是否成立，并据 P1-5 决定 Agent 是否升级
- [ ] 论文初稿与可复现包

---

## 7. 评估方案

### 设计空间是二维的

揭示策略与最终发布策略必须分开，二者正交。默认"只公开最终 p 值即为低泄露"是一个需要被检验的假设，而不是设计前提。

维度一：最终发布策略（final release）

| 档位 | 公开内容 | 在实验中的角色 | 预期泄露 |
|---|---|---|---|
| F1 | 仅通过显著性阈值的命中 | 最保守的现实发布 | 最低 |
| F2 | top-k 命中位点 | 常见论文发布形式 | 低 |
| F3 | 全部位点 p 值 | 逐步放宽 | 中 |
| F4 | 全部位点 p 值 + 效应量 + 标准误 | 接近完整汇总统计 | 高 |
| F5 | 完整汇总统计 | 泄露上界与攻击 sanity check | 最高 |

维度二：中间揭示策略（intermediate reveal）

| 档位 | 公开内容 | 在实验中的角色 | 预期边际泄露 |
|---|---|---|---|
| M0 | 无诊断量 | 下界 | — |
| M1 | 仅 λ_GC | 最常见的实践 | 待测 |
| M2 | 主成分特征值谱与主成分-表型相关 | 分层诊断所需 | 待测 |
| M3 | 各位点缺失率与 QC 统计 | 流程质控所需 | 待测 |

P0-4 在固定 F 档的前提下扫描 M 档，测的是 M 相对 M0 的增量，且至少在保守档与宽松档各测一次；P1-4 做全面交叉。

### 指标定义

公开对象（句法层）。不使用单一比特数或条目数。每条公开记录按定义 1 的八元组落盘。报告四个量：承诺集合大小、单次实际公开集合大小、跨运行账本规模、以及账本上被标记的关系数。句法量不对应任何隐私声明，只用于审计与可复现。

一致性（句法层）。允许子集一致性与精确路径一致性分别报告，前者是安全性质，后者是完整性性质，两者不可互换（定义 2）。

泄露（经验层）。两个指标，均需报告对手先验档位：

- 成员推断：报告 TPR@FPR≤0.1% 与 TPR@FPR≤1% 两档，不使用平均准确率。至少两类攻击族分别报告（聚合等位基因频率类与似然比/ML 类）。必须逐项写清：member/non-member 的构造方式、两者是否祖先匹配、目标个体是否出现在攻击者持有的参考面板中、面板重叠与不重叠两档、攻击者是否知道病例对照状态、以及逐祖先分别的 TPR。任何 TPR 提升都必须在祖先匹配且面板不重叠的严格控制下复现，否则不得归因于发布策略。
- 基因型恢复：全文只使用一个主定义，建议取"恢复单个目标个体在单个位点上的基因型类别（0/1/2）"，其余口径（等位基因存在性、minor allele count）只在附录作敏感性分析，绝不混用。报告可恢复的变异比例。

精度格（precision lattice）。四元组中的 precision 不是一个标量，需要按量的类型分别定义：对 p 值定义有效数字位数与分箱边界；对效应量与标准误定义有效数字位数与舍入方式；对等位基因频率定义分箱；对缺失率与 QC 统计定义分箱；对诊断量定义 top-k 截断与阈值布尔。跨量的精度之间只定义偏序（更粗 / 不可比），不假设存在统一刻度——Wagner 的非单调性结论正是禁止这种假设的依据。

账本关系的形式化。定义 3 的四类关系落在同一 quantity 上，按如下方式判定：细化定义为精度格上严格更粗→更精的序关系；条件放宽定义为 condition 从 conditional(·) 变为 always；重复定义为 quantity 与 precision 均相同；接收方扩张定义为接收方集合严格超集。并集覆盖定义为：若若干次发布的并集，对某个 quantity 所能提供的信息不少于某次单次更高精度发布的至少一个可见档位，则触发标记。所有判定都只在句法层进行，不对攻击风险作任何声明。

效用。功效与假阳性率需明确区分模拟真信号位点与 null 位点，并对位点做 LD 剪枝或按 LD block 处理，否则位点间相关性会污染功效与 FPR 的解释，不能把 LD 相关位点当作独立样本逐 SNP 计算 FPR。λ_GC 作为校准指标报告，但需说明其局限，不单独用它判定分层已被校正。

强制边界。报告三张清单：TCB 覆盖范围内被阻止的违约、被检测但未阻止的违约、落在 TCB 之外的已知行为；同时报告 TCB 规模（组件数、代码行数）与 logging-only 基线在同样注入下的表现。

开销。报告 MPC 墙钟时间与通信量随位点数、站点数的伸缩曲线，并写明硬件与网络条件。不追求绝对速度。

编译质量。报告安全关键字段错误率（漏报揭示、多报揭示、精度错误、接收方错误、条件错误、路径公开量低估），不使用 JSON 完全匹配。

### 基线方法

1. 合池明文金标准（PLINK2）— 效用上界，非隐私基线。
2. 联邦 meta 分析四档（固定效应 / 随机效应 / 分层 / 含 genomic control）— 现实工程主流，避免 strawman。
3. 完整汇总统计发布 — 泄露上界，攻击实现的 sanity check。
4. Pascoal / GenDPR 风格的安全发布风险评估 — 最重要的领域相关基线。不要求完整复现其系统，但必须对齐其攻击设定与安全发布条件，并逐项说明本工作覆盖了它不能覆盖的哪些场景（运行期中间揭示与跨运行累积）。
5. logging-only 工作流日志与 Nextflow/CWL 风格插桩的联邦 pipeline — 回答"与普通 provenance 系统有何不同"。
6. 无预承诺的自适应流程 — 用于量化预承诺的统计代价。
7. 规则/模板编译器 — Agent 必须击败的强基线，仅与"固定脚本"比较不够。
8. 无强制机制的事后审计 — 只比对日志、不拦截揭示，用于隔离强制机制本身的贡献。
9. 受控分析平台（DataSHIELD 类）与 DP 查询平台（OpenDP 类）— 前者作为"限制操作集合"路线的对照，后者作为"我们不提供 ε-accounting"的对照；两者都只做能力对照，不实现。
10. TEE 路线 — 仅做信任假设对照，不实现。

### 统计分析

- 主实验 3 次独立重复，报告均值 ± 标准差。
- 功效与 FPR 用 bootstrap（1000 次重采样）给出置信区间。
- 编译质量对比用 McNemar 配对检验。
- 二维扫描中的趋势用 Jonckheere–Terpstra 有序趋势检验。
- 所有结论报告效应量，不只报显著性。

---

## 8. 预期输出

```
gwas-reveal-governance/
├── code/{data,baseline,commit,runtime,mpc,audit,attacks}/
├── bench/compilebench/     # 意图→揭示策略 基准集（默认附录）
├── results/{figures,tables,raw_data}/
└── audit_reports/          # 审计报告样例
```

对外交付三件：开源仓库（3 站点模拟数据上一条命令跑通）；TCB 清单与违约注入测试集；论文。

论文承诺 C1–C4，重心在 C3。Agent 前端与编译质量基准默认位于附录，仅当 P1-5 显示显著优势时才考虑升级。

---

## 9. 风险与备选方案

- **风险 1｜最终发布主导泄露，中间揭示的边际效应可忽略。** 这是本工作最可能被推翻的假设，且 P0-4 会最先给出答案。
  → 备选：论文改写为**多站点 GWAS 工作流的发布策略经验前沿**——比较 F1–F5 的隐私风险与效用，中间揭示作为负结果保留。动机从"中间揭示的预算"变为"发布策略的选择本身缺乏可审计依据"。这条路依然完整，且负结果本身有引用价值。
- **风险 2｜系统机制被判定为"领域化的 declassification allowlist + provenance log"。** 这是新颖性层面的主要威胁。
  → 备选：把论文重心从机制彻底转向测量——强调这是第一个把揭示策略、强制边界、跨运行审计与基因组攻击评估连成可复现系统的多站点 GWAS 工作；交付 TCB 清单、违约注入测试集、logging-only 与联邦插桩两个对照基线、以及公开策略前沿数据。同时用 P0-5 的对照基线与 P0-7 的账本检测，正面回答"和普通审计日志有何不同"。标题避免"新框架"一类措辞，改用"强制原型与经验研究"的定位。
- **风险 3｜公开集合与实际攻击风险不相关（度量非单调）。** Wagner 的工作表明这在基因组隐私度量中是常见现象，必须视为真实风险而非小概率。
  → 备选：放弃任何单一标量的"预算"叙事，改为攻击模型特定的风险报告。论文标题与主张从 privacy budget 改为 empirical leakage reporting，这与本工作已有的双层立场一致，属于顺水推舟的降级。

- **风险 4｜强制机制挡不住现实的绕过，A2 的覆盖范围被实测证伪。**
  → 备选：把 A2 进一步收窄为"只能提交由清单生成的 DSL 工件、不能提交任意代码"，把任意代码环境列为 out-of-scope 与 future work；C1 相应地从"执行不变式"降级为"生成工件上的 trace 合规性质"。用 P0-5 把这个边界测出来并写进论文，比含糊声称全覆盖要好得多。
- **风险 5｜Agent 在 CompileBench pilot 上不优于模板基线。**
  → 备选：Agent 移入附录或降为可选接口。主贡献保持 C1–C4 不变，论文反而更容易被安全社区接受。这条路线已在 1.5 中事前声明。
- **风险 6｜预承诺的统计代价过大。**
  → 备选：把负结果转为贡献——量化多站点 GWAS 中预承诺分析的统计代价，并讨论何种分析适合预承诺、何种需要噪声化诊断量。需提前与导师对齐，避免把负结果误判为实验失败。
- **风险 7｜攻击实现被质疑太弱，或攻击增益被祖先结构混淆。** 二者都会让"泄露降低"的结论不可信。
  → 对策：把 F5 端点、双攻击族、双先验，以及祖先匹配、参考面板重叠两档、逐祖先 TPR 全部作为攻击强度的验收门槛（P0-3 的 Pass 条件）。任一项不过，不允许进入任何隐私结论，只能报告"在这些攻击下未观测到可测增加"。
- **风险 8｜工作量过大，四条线中至少一条只有半成品。** 明文基线、攻击、强制机制、MPC 绑定任何一条薄弱都会削弱可信度；安全主成分尤其容易在 MPC 中过慢或数值不稳。
  → 对策：以 P0 的八项作为最小可发表集，其余（预承诺代价、扩展版账本与消融、二维扫描、CompileBench）在 P0 通过后再推进。P0-6 从设计上就按"安全聚合等位基因计数 + 逐位点关联统计量"的最小形态做，不做完整安全 PCA；超预算即触发降级为明文模拟揭示策略加小规模绑定演示。
- **风险 9｜导师期待交付可部署系统，而本工作交付的是带基准的研究。**
  → 早沟通。若目标确为部署，需重新调整 C1–C4 与附录项的权重，这应在开题前当面确认。

---

## 10. 工作量估计

- 环境准备 · 中等 — 难点是补 GWAS 领域知识，而非装环境。以"用 PLINK2 亲手跑通一次完整流程"为硬性验收。
- 核心实现 · 重度 — 关键路径是强制机制（DSL、中介运行时、日志链）、执行层与攻击实现。明文基线与承诺清单 schema 必须先行且做扎实，它们是后续一切的参照系。
- 实验运行 · 中等偏重 — 瓶颈是 MPC 重复运行时间与二维扫描的组合爆炸，必须提前规划机时并坚持 scope 纪律。
- 结果分析 · 中等 — 三张核心图加三张表；难点在负结果的解释与主张降级的决策时机。
- 总计 · ⭐⭐⭐⭐ — 完整方案约 6–9 个月全职。若只保留最小可行版本（明文金标准、四档联邦基线、双攻击族、强制边界与违约注入、发布策略前沿、账本 P0 部分、小规模 MPC demo），约 4–5 个月。

里程碑

- M1：明文全流程与站点切分数据就绪（P0-1）。
- M2：联邦 meta 四档基线的失效实证出结果（P0-2）。这个数字决定全篇说服力。
- M3：攻击端点在 F5 上通过验收（P0-3）。
- M4：中间揭示边际泄露出结果（P0-4）。此处决定论文走哪条路。
- M5：承诺清单、强制机制、TCB 清单与绕过注入完成，含两个对照基线（P0-5）。此处的对照结果决定"MPC 更适合"这一假设是否保留。
- M6：MPC 揭示绑定演示与金标准对齐，或触发降级（P0-6）。
- M7：跨运行释放关系与策略字段消融的最小版出结果（P0-7、P0-8）。
- M8：预承诺代价、扩展版账本与消融、二维扫描、CompileBench pilot 完成，三张核心图成型（P1-1 至 P1-5）。

---

## 附：开题与答辩需当场回答的七个问题

1. 为什么用 MPC 作为执行底座？这是一个待检验假设，不是既定前提：MPC 的 reveal / open 是数量有限、位置明确的显式原语，满足"唯一外部输出且 I/O 被中介"的条件时，就可以成为窄的强制点。我们把它写成假设，并用 P0-5 的两个对照基线（logging-only provenance 与联邦插桩 pipeline）在同一组绕过注入下实测。若三条线表现接近，本工作改为"MPC 是一个清晰的实例化后端"，贡献收回到领域语义清单与经验评估。
2. 为什么不是 TEE？信任假设不同。TEE 把信任转移给硬件厂商与 enclave 实现；MPC 不需要该假设。只做定性对照。
3. 你的系统到底防谁？防受框架约束、无法执行任意代码、无法绕过中介运行时的作者（A2）；不防能改写执行框架或篡改日志的作者。这条边界由 TCB 清单显式界定，并由 P0-5 逐条实测。
4. 这和一份 declassification allowlist 加一个日志有什么区别？区别有三处，且都可测：揭示边界是被强制的唯一出口，而不是被记录的副作用；判定绑定在揭示对象语义上（精度、接收方、条件、推导来源）；风险估计绑定在具体基因组攻击模型上。1.6 的能力对照表逐格列出差异，P0-5 的 logging-only 与联邦插桩基线、P1-1 的账本检测分别把这些差别量化出来。
5. 为什么需要 Agent？用模板编译器不行吗？这是待验证问题，不是前提，且 Agent 不在标题与主贡献中。P1-5 用安全关键错误率对比，若不显著优于模板基线，Agent 停留在附录。
6. Agent 会不会自己成为泄露点？Agent 位于数据平面之外，其输入不含任何私有统计量，输出只经过 schema 校验而不直接触达执行层。但这只说明它不直接观察私有输入，不说明它生成的策略是好的——后者只能靠经验层评测。
7. 数据从哪来？1000 Genomes 公开人类基因组资源与 msprime 模拟数据；不使用受控访问的临床队列，也不新收集受试者数据。应主动说明，这是合规与伦理上的实质优势。
