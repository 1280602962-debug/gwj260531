# SECONDARY SCORING SOFTWARE AUDIT

Protocol: SCORING_AND_DOCKING_METHOD_ABLATION_V2

- GNINA binary: `/mnt/d/CADD paper exercise/gnina/bin/gnina`
- GNINA version: `gnina v1.3.2 master:f23dd2b   Built Jul  8 2025.`
- activate: `/mnt/d/CADD paper exercise/gnina/activate.sh` exists=True
- libcudnn: `True`
- help cnn_scoring default rescore: True
- help pose_sort_order CNNscore default: True
- CNN model: built-in default ensemble when `--cnn` / `--cnn_model` omitted (frozen)
- frozen flags: `--no_gpu --scoring vina --cnn_scoring rescore --pose_sort_order CNNscore`
- forbidden: `--minimize`, `--cnn_scoring all`

- RTMScore python: `/home/gwj/miniconda3/envs/rtmscore/bin/python`
- python/torch:
```
3.10.20
2.6.0
```
- checkpoint: `/home/gwj/software/RTMScore/trained_models/rtmscore_model1.pth` sha256=399727a85de96c011a75a4c242ebe947753724f01e8d15139fae7bcad106bfc2
- model2/model3: **not loaded, not compared**
- load model1: rc=0 out=`dict 2` err=``
- example infer rc1=0 rc2=0
- example numeric gate: none (identity only; no required score value)

- PHASE1_PASS: **YES**
