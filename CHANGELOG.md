# Changelog

## v1.3.0
- macOS/Homebrew 环境自动探测 `hmmscan`。
- 自动把 HMMER bin 目录注入当前 PyMOL 进程 PATH。
- `cdr_backend` 显示实际 hmmscan、Python 和识别后端信息。

## v1.2.0
- VHH / 单 VH 仅创建 H_CDR1/2/3 三个 selection。
- 常规 VH+VL 创建 H/L 各三个 selection。
- 删除旧版重复的对象级汇总 selection。
- 重新识别前自动清理旧版 `abCDR_*` selection。

## v1.1.0
- 中文 GUI。
- CDR1/2/3 增加颜色下拉选择。
- 改进 VHH heuristic CDR3 识别，允许 CDR3 内额外 Cys。

## v1.0.0
- 初始版本。
