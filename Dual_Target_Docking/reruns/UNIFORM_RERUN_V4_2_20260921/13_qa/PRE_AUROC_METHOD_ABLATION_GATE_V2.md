# PRE_AUROC_METHOD_ABLATION_GATE_V2

PRE_AUROC_GATE = PASS

## Answers
### 1. 当前脚本的配置权威是否唯一？
YES — STALE_PROVENANCE_ONLY; runners use ablation_config.py  (internal_ok=True)

### 2. physical / pair-expanded计数是否区分正确？
YES — pair-expanded 1614; physical 1592 / 14250  (internal_ok=True)

### 3. M0是否锁定？
YES  (internal_ok=True)

### 4. M2是否锁定？
YES  (internal_ok=True)

### 5. 原61个M1失败是否全部分类？
YES — 61 GNINA_INPUT_TYPE_FAILURE, n_retries=1  (internal_ok=True)

### 6. M1是否已采用统一、保持真实化学图和原始Vina坐标的GNINA-compatible representation？
YES  (internal_ok=True)

### 7. M1 physical pose rows是否应为14250，且实际全部入账或明确缺失？
long=14097  (internal_ok=True)

### 8. M1是否严格来自Vina mode1？
YES  (internal_ok=True)

### 9. M1b是否只从同一Vina saved pose set取max？
YES  (internal_ok=True)

### 10. M3 parser bug是否已经从原始out.pdbqt修复？
YES — MASTER_REPARSED  (internal_ok=True)

### 11. 所有现有合法M3 SUCCESS是否具有finite mode1 CNNscore？
YES  (internal_ok=True)

### 12. M3宏环/特殊atom representation是否完成全量审计？
YES — 29 ligands, CGn≠Gn  (internal_ok=True)

### 13. 若存在无效M3输入，是否只按预先冻结affected set统一修复、没有按结果择优？
YES  (internal_ok=True)

### 14. AB_040是否未作为独立observation？
YES  (internal_ok=True)

### 15. shared receptor-parent jobs是否正确跨pair复用？
YES — 22 shared physical jobs  (internal_ok=True)

### 16. activity_eligible=0是否不会进入four-class？
YES — class empty, not neither  (internal_ok=True)

### 17. 是否存在任何score imputation？
NO  (internal_ok=True)

### 18. 是否存在任何post hoc model/seed/receptor/box/pose/engine selection？
NO  (internal_ok=True)

### 19. 是否修改了任何原始文件？
NO  (internal_ok=True)

### 20. 本任务是否计算/查看了正式AUROC、summary_min、ΔAUC、Top-K、EF？
NO  (internal_ok=True)

No AUROC was computed.
