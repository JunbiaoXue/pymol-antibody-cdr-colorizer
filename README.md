# PyMOL Antibody CDR Colorizer

一个支持 Python 3 的 PyMOL 抗体 CDR 自动识别与着色插件，支持常规 VH/VL 和 VHH/Nanobody。

> 当前版本：v1.3.0 · 识别优先级：ANARCI → ANARCII → 内置 heuristic · CDR 定义：IMGT


## v1.3.0：macOS / Homebrew 的 hmmscan 自动修复

macOS 从 Finder/Launchpad 启动 `PyMOL.app` 时，GUI 进程通常不会继承终端的 Homebrew `PATH`。因此即使终端执行 `which hmmscan` 能看到 `/opt/homebrew/bin/hmmscan`，ANARCI 仍可能报：

```text
[Errno 2] No such file or directory: 'hmmscan'
```

v1.3.0 会在插件启动及调用 ANARCI 前自动搜索并注入：

```text
当前 PATH
PyMOL Python / sys.prefix 的 bin
/opt/homebrew/bin/hmmscan
/usr/local/bin/hmmscan
/opt/local/bin/hmmscan
Conda / Mamba / virtualenv 的 bin
```

只修改**当前 PyMOL 进程**的 `PATH`，不会修改 `.zshrc`、系统环境或 PyMOL.app 文件。

可在 PyMOL 中运行：

```text
cdr_backend
```

查看插件实际使用的后端、`hmmscan` 绝对路径和 PyMOL Python 路径。

## 功能

- 打开 PDB/mmCIF 后自动扫描抗体可变域。
- 支持常规抗体 **VH/VL**，同时支持 **VHH / Nanobody 单域抗体**。
- 按 **IMGT** 定义显示 CDR：
  - CDR1：27–38
  - CDR2：56–65
  - CDR3：105–117
- CDR1/CDR2/CDR3 分别着色。
- 中文 GUI。
- 每个 CDR 的颜色使用下拉选项框选择。
- 可选将 CDR 显示为 sticks。
- 自动建立 CDR selections，方便后续 PyMOL 操作。

## VHH 兼容性

v1.2.0 专门修正了 VHH fallback 识别。

旧逻辑在没有 ANARCI/ANARCII 时，会把 FR4 `WGQG/FGxG` 前的**最后一个 Cys**当成 CDR3 起点。这对普通 VH 经常可用，但对 VHH 不可靠，因为 VHH 的 CDR3 较长，并且经常包含额外 Cys/非经典二硫键。

新版改为：

1. 先定位 C 端 FR4 motif；
2. 再寻找 FR3 末端的保守 Cys（IMGT 104 附近）；
3. 允许 CDR3 内继续出现一个或多个 Cys；
4. 不要求同时存在轻链，因此单独加载 VHH 也能识别；
5. 不再强依赖 PDB 的 H/L 链名，优先根据序列判断重链/轻链。

> 注意：内置 heuristic 仍然是近似识别。需要正式 IMGT 编号时，建议在 PyMOL 的 Python 环境中安装 ANARCI 或 ANARCII。

## 安装

PyMOL 中：

1. `Plugin` → `Plugin Manager`
2. `Install New Plugin`
3. 选择插件 ZIP 文件
4. 安装完成后重新打开 PyMOL（建议）
5. 菜单中选择：`Plugin` → `抗体 CDR 自动识别与着色`

## GUI

插件面板包含：

- 结构对象选择
- CDR1 颜色下拉框
- CDR2 颜色下拉框
- CDR3 颜色下拉框
- 自动处理新加载结构
- CDR sticks 显示
- “识别并着色”
- “扫描全部结构”
- “删除 CDR 选择”

默认颜色：

- CDR1：黄色 `yellow`
- CDR2：洋红 `magenta`
- CDR3：青色 `cyan`

## PyMOL 命令

仍保留命令行接口：

```text
cdr_color
cdr_color my_object
cdr_color all, red, green, blue, 1
cdr_info my_object
cdr_auto 1
cdr_auto 0
cdr_backend
cdr_clear my_object
```

其中：

```text
cdr_color selection, CDR1颜色, CDR2颜色, CDR3颜色, sticks
```

例如：

```text
cdr_color my_vhh, yellow, magenta, cyan, 1
```

## 识别后端

优先级：

```text
ANARCI → ANARCII → 内置 heuristic
```

运行：

```text
cdr_backend
```

可以查看当前 PyMOL 实际使用的识别后端和 Python 路径。

### ANARCI / ANARCII

如果安装了编号工具，插件会基于结构提取出的氨基酸序列做 IMGT 编号，不要求结构链名一定是 `H` / `L`。

VHH 会作为重链型可变域处理，因此**不需要配对轻链**。

## 插件生成的 selection

v1.2.0 起取消重复的“汇总 CDR”selection，只保留实际抗体链对应的 CDR。

单链 VHH / 单独 VH（只检测到一个 H 型可变域）只生成 3 个：

```text
abCDR_myVHH_H_CDR1
abCDR_myVHH_H_CDR2
abCDR_myVHH_H_CDR3
```

常规 VH + VL 抗体生成 6 个：

```text
abCDR_myFab_H_CDR1
abCDR_myFab_H_CDR2
abCDR_myFab_H_CDR3
abCDR_myFab_L_CDR1
abCDR_myFab_L_CDR2
abCDR_myFab_L_CDR3
```

其中 κ、λ 都统一显示为 `L`。如果同一对象中存在多个 H 或多个 L 可变域，插件会自动把原始 chain ID 加入 selection 名称，避免重名。

每次重新识别时，会先删除该对象由旧版插件创建的 `abCDR_*` selection，因此从 v1.0/v1.1 升级后不会残留旧的重复 selection。

## 局限

- heuristic 模式不能替代正式抗体编号工具。
- 极端工程化 scaffold、严重缺失残基或结构中只保留局部片段时可能无法自动识别。
- PyMOL 结构中的残基编号不一定等于 IMGT 编号；插件是先基于序列识别 CDR，再映射回 PyMOL 结构残基。
