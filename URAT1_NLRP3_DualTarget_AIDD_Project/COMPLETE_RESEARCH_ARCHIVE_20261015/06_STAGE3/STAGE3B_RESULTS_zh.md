# 第三阶段 B

原始报告仍在 `results/nlrp3_model_v2/stage3b_local/analysis/STAGE3B_REPORT_zh.md`。这里只保留已经和账本、汇总表核对过的数字。

194 次任务全部 `completed`，每次 9 个姿势，没有失败，没有超时。单次 67.4–280.8 秒，中位数 146.3 秒。并行是 2 个任务乘 4 个线程。软件、盒子、穷举度和评分与生产协议相同，受体文件与历史生产受体字节相同。

32 个配体全部准备成功。形式电荷 −1 的 24 个，−2 的 3 个，0 的 5 个。6 个有未指定立体中心，没有补指派。

96 个 URAT1 评价和 96 个 NLRP3 评价全部对接完成。URAT1 A2 通过 88 个种子，8 个种子没有几何合格姿势。NLRP3 结构门通过 65 个种子。同种子双通过 62 个。

分子层面：13 个结构合格且化学标记干净，8 个结构合格但化学标记不通过，11 个双通过种子少于 2 个。13 个是 ARONIXIL、BENFURODIL HEMISUCCINATE、BROXITALAMIC ACID、CEDIROGANT、GSK-065、GSK-625433、INDEGLITAZAR、MEBROFENIN、METANIXIN、PIRINIXIC ACID、RECOFLAVONE、TETOMILAST（各 3 个双通过种子）和 FEBUXOSTAT（2 个）。8 个里 RAZUPENEM、FENBENICILLIN、TEBIPENEM 是 β-内酰胺；SUXIBUZONE、OLSALAZINE、N6022、ENOXAMAST、SUCCISULFONE 的 PAINS/Brenk 标记不干净。11 个里多数是 URAT1 已通过而 NLRP3 未通过，不是对接失败，也不是“没有 URAT1 活性”。GSK-2636771 是唯一一个三个种子上两个靶点都没有合格姿势的分子，它仍然写出了 9 个姿势。

REP_07837（MK-8457）的两次补算各 9 个姿势。URAT1 三个种子的 A2 通过。NLRP3 三个种子的最高 CNNscore 约 0.27–0.28，结构门未过。双通过种子为 0。历史日志里的超时记录没有被改写成失败的几何评价。

已存 156 个分子重读后，冻结的 33 个候选全部仍然结构合格且化学标记通过。另有 21 个已存分子达到结构门，但不在冻结化学表里，没有自动提名。新的 13 个与这 33 个没有 ID 重叠。

PF-03882845 三个种子都双通过。URAT1 羧酸氧到 Arg477 约 3.23–3.24 Å，质心约 1.33 Å，CNNscore 约 0.643。NLRP3 重叠 0.867，IFP 0.85，7/7 关键接触，无冲突。没有用分子动力学改这些判定。
