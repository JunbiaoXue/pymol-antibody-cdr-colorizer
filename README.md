# PyMOL Antibody CDR Colorizer

一个支持 Python 3 的 PyMOL 抗体 CDR 自动识别与着色插件，支持常规 VH/VL 和 VHH/Nanobody。

> 当前版本：v1.3.0 · 识别优先级：ANARCI → ANARCII → 内置 heuristic · CDR 定义：IMGT

## 功能

- 打开 PDB/mmCIF 后自动扫描抗体可变域。
- 支持常规抗体 **VH/VL**，同时支持 **VHH / Nanobody 单域抗体**。
- 按 **IMGT** 定义显示 CDR：CDR1 27–38、CDR2 56–65、CDR3 105–117。
- CDR1/CDR2/CDR3 分别着色，中文 GUI，颜色下拉选择。
- VHH / 单独 VH 仅创建 3 个 H_CDR selection；VH+VL 创建 H/L 各 3 个。
- 优先使用 ANARCI；无编号工具时提供内置 heuristic fallback。
- macOS v1.3.0 自动探测 Homebrew HMMER `hmmscan` 并注入 PyMOL 当前进程 PATH。

完整说明、源码与可安装 ZIP 将随 v1.3.0 发布提交加入本仓库。

## License

MIT License © 2026 JunbiaoXue
