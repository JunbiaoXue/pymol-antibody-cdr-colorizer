"""PyMOL 抗体 CDR 自动识别与着色插件。"""
from __future__ import annotations

__version__ = "1.3.0"
__author__ = "JunbiaoXue"

_DIALOG = None
_WATCHER = None
_AUTO_ENABLED = True
_LAST_SIGNATURES = {}
_LAST_RESULTS = {}

# 中文显示名 -> PyMOL 内置颜色名
_COLOR_OPTIONS = [
    ("黄色", "yellow"),
    ("洋红", "magenta"),
    ("青色", "cyan"),
    ("红色", "red"),
    ("绿色", "green"),
    ("蓝色", "blue"),
    ("橙色", "orange"),
    ("紫色", "violet"),
    ("鲑红", "salmon"),
    ("柠檬绿", "lime"),
    ("海蓝", "marine"),
    ("蓝绿色", "teal"),
    ("天蓝", "skyblue"),
    ("白色", "white"),
    ("浅灰", "gray70"),
]


def _prefs_get(key, default):
    try:
        from pymol.plugins import pref_get
        return pref_get("antibody_cdr_colorizer_" + key, default)
    except Exception:
        return default


def _prefs_set(key, value):
    try:
        from pymol.plugins import pref_set, pref_save
        pref_set("antibody_cdr_colorizer_" + key, value)
        pref_save()
    except Exception:
        pass


def _object_signature(name):
    from pymol import cmd
    try:
        return (cmd.count_atoms(f"%{name}"), cmd.count_states(f"%{name}"))
    except Exception:
        return None


def _eligible_objects():
    from pymol import cmd
    names = []
    for name in cmd.get_names("objects"):
        try:
            if cmd.get_type(name) == "object:molecule" and cmd.count_atoms(f"%{name} and polymer.protein") > 0:
                names.append(name)
        except Exception:
            continue
    return names


def _backend_cn(backend):
    if backend == "heuristic":
        return "内置启发式"
    return backend


def _color_one(object_name, quiet=False, show_sticks=None):
    from .core import analyze_object, apply_colors, format_result

    if show_sticks is None:
        show_sticks = bool(_prefs_get("show_sticks", False))
    result = analyze_object(object_name)
    _LAST_RESULTS[object_name] = result

    if result.domains:
        apply_colors(
            result,
            cdr1_color=str(_prefs_get("cdr1_color", "yellow")),
            cdr2_color=str(_prefs_get("cdr2_color", "magenta")),
            cdr3_color=str(_prefs_get("cdr3_color", "cyan")),
            show_sticks=bool(show_sticks),
        )
    if not quiet:
        print("[抗体 CDR 着色]\n" + format_result(result))
        if result.backend == "heuristic":
            print("[抗体 CDR 着色] 当前使用内置启发式识别；正式 IMGT 编号建议安装 ANARCI/ANARCII。")
    return result


def cdr_color(selection="all", cdr1="yellow", cdr2="magenta", cdr3="cyan", sticks=0, quiet=0):
    """自动识别并给抗体 CDR1/2/3 着色。

    PyMOL 命令示例：
        cdr_color
        cdr_color my_object
        cdr_color all, yellow, magenta, cyan, 1
    """
    from pymol import cmd

    _prefs_set("cdr1_color", str(cdr1))
    _prefs_set("cdr2_color", str(cdr2))
    _prefs_set("cdr3_color", str(cdr3))
    _prefs_set("show_sticks", bool(int(sticks)))

    objects = cmd.get_object_list(f"({selection})") if selection not in ("all", "*") else _eligible_objects()
    count = 0
    for name in objects:
        result = _color_one(name, quiet=bool(int(quiet)), show_sticks=bool(int(sticks)))
        if result.domains:
            count += 1
        _LAST_SIGNATURES[name] = _object_signature(name)
    return count


def cdr_info(selection="all"):
    """打印检测到的抗体可变域、CDR 序列和结构残基编号。"""
    from pymol import cmd
    from .core import format_result

    objects = cmd.get_object_list(f"({selection})") if selection not in ("all", "*") else _eligible_objects()
    for name in objects:
        result = _LAST_RESULTS.get(name) or _color_one(name, quiet=True)
        print("[抗体 CDR 着色]\n" + format_result(result))


def cdr_clear(selection="all"):
    """删除插件创建的 CDR selections；不恢复已有原子颜色。"""
    from pymol import cmd
    from .core import delete_plugin_selections

    if selection in ("all", "*"):
        delete_plugin_selections(None)
    else:
        for name in cmd.get_object_list(f"({selection})"):
            delete_plugin_selections(name)


def cdr_auto(on=1):
    """开启/关闭新加载结构的自动识别和着色。"""
    global _AUTO_ENABLED
    _AUTO_ENABLED = bool(int(on))
    _prefs_set("auto_enabled", _AUTO_ENABLED)
    print(f"[抗体 CDR 着色] 自动识别：{'开启' if _AUTO_ENABLED else '关闭'}")
    return int(_AUTO_ENABLED)


def cdr_backend():
    """显示当前编号后端、hmmscan 路径及 PyMOL Python 环境。"""
    import os
    import sys
    from .core import backend_status, hmmscan_status

    hmmscan = hmmscan_status()
    backend = backend_status()
    print(f"[抗体 CDR 着色] 当前后端：{_backend_cn(backend)}")
    print(f"[抗体 CDR 着色] hmmscan：{hmmscan or '未找到'}")
    print(f"[抗体 CDR 着色] PyMOL Python：{sys.version.split()[0]}  路径：{sys.executable}")
    if hmmscan:
        print(f"[抗体 CDR 着色] HMMER bin 已加入当前 PyMOL PATH：{os.path.dirname(hmmscan)}")
    elif backend != "ANARCII":
        print("[抗体 CDR 着色] 未找到 hmmscan；ANARCI 无法运行。请确认 HMMER 已安装。")
    if backend == "heuristic":
        print("[抗体 CDR 着色] 当前将使用内置启发式识别。")
    return backend


def _watch_tick():
    if not _AUTO_ENABLED:
        return
    current = set(_eligible_objects())

    for old in list(_LAST_SIGNATURES):
        if old not in current:
            _LAST_SIGNATURES.pop(old, None)
            _LAST_RESULTS.pop(old, None)

    for name in current:
        sig = _object_signature(name)
        if sig is None:
            continue
        if _LAST_SIGNATURES.get(name) == sig:
            continue
        try:
            _color_one(name, quiet=True)
        except Exception as exc:
            print(f"[抗体 CDR 着色] 自动扫描 {name} 失败：{exc}")
        finally:
            _LAST_SIGNATURES[name] = sig


def _start_watcher():
    global _WATCHER, _AUTO_ENABLED
    from pymol.Qt import QtCore

    _AUTO_ENABLED = bool(_prefs_get("auto_enabled", True))
    if _WATCHER is None:
        _WATCHER = QtCore.QTimer()
        _WATCHER.setInterval(1200)
        _WATCHER.timeout.connect(_watch_tick)
        _WATCHER.start()
    QtCore.QTimer.singleShot(300, _watch_tick)


def _fill_color_combo(combo, current_color):
    combo.clear()
    for cn_name, pymol_name in _COLOR_OPTIONS:
        combo.addItem(f"{cn_name}  ({pymol_name})", pymol_name)
    idx = combo.findData(str(current_color))
    if idx < 0:
        # 保留旧版本或用户脚本保存的自定义 PyMOL color 名。
        combo.addItem(f"自定义  ({current_color})", str(current_color))
        idx = combo.count() - 1
    combo.setCurrentIndex(idx)


def run_plugin_gui():
    global _DIALOG
    from pymol.Qt import QtWidgets, QtCore
    from .core import backend_status, format_result, delete_plugin_selections

    if _DIALOG is not None:
        _DIALOG.show()
        _DIALOG.raise_()
        _DIALOG.activateWindow()
        return

    dialog = QtWidgets.QDialog()
    dialog.setWindowTitle("抗体 CDR 自动识别与着色")
    dialog.resize(620, 500)
    layout = QtWidgets.QVBoxLayout(dialog)

    intro = QtWidgets.QLabel(
        "自动识别普通抗体 VH/VL 及 VHH/纳米抗体可变域，并按 IMGT CDR1/2/3 分别着色。\n"
        "VHH/单独 VH 仅创建 3 个 H_CDR selection；VH+VL 创建 6 个 H/L CDR selection。"
    )
    intro.setWordWrap(True)
    layout.addWidget(intro)

    backend_label = QtWidgets.QLabel(f"当前编号后端：{_backend_cn(backend_status())}")
    layout.addWidget(backend_label)

    form = QtWidgets.QFormLayout()
    object_box = QtWidgets.QComboBox()
    form.addRow("结构对象：", object_box)

    c1 = QtWidgets.QComboBox()
    c2 = QtWidgets.QComboBox()
    c3 = QtWidgets.QComboBox()
    _fill_color_combo(c1, _prefs_get("cdr1_color", "yellow"))
    _fill_color_combo(c2, _prefs_get("cdr2_color", "magenta"))
    _fill_color_combo(c3, _prefs_get("cdr3_color", "cyan"))
    form.addRow("CDR1 颜色：", c1)
    form.addRow("CDR2 颜色：", c2)
    form.addRow("CDR3 颜色：", c3)
    layout.addLayout(form)

    opts = QtWidgets.QHBoxLayout()
    auto_cb = QtWidgets.QCheckBox("自动处理新加载的结构")
    auto_cb.setChecked(_AUTO_ENABLED)
    sticks_cb = QtWidgets.QCheckBox("CDR 显示为 sticks")
    sticks_cb.setChecked(bool(_prefs_get("show_sticks", False)))
    opts.addWidget(auto_cb)
    opts.addWidget(sticks_cb)
    opts.addStretch(1)
    layout.addLayout(opts)

    buttons = QtWidgets.QHBoxLayout()
    apply_btn = QtWidgets.QPushButton("识别并着色")
    scan_btn = QtWidgets.QPushButton("扫描全部结构")
    clear_btn = QtWidgets.QPushButton("删除 CDR 选择")
    buttons.addWidget(apply_btn)
    buttons.addWidget(scan_btn)
    buttons.addWidget(clear_btn)
    layout.addLayout(buttons)

    status = QtWidgets.QPlainTextEdit()
    status.setReadOnly(True)
    status.setPlaceholderText("识别结果会显示在这里。")
    layout.addWidget(status, 1)

    note = QtWidgets.QLabel(
        "IMGT 定义：CDR1 27–38，CDR2 56–65，CDR3 105–117。"
        "无 ANARCI/ANARCII 时会自动使用内置启发式模式；新版已针对 VHH 长 CDR3 和额外 Cys 做兼容。"
    )
    note.setWordWrap(True)
    layout.addWidget(note)

    def refresh_objects():
        current = object_box.currentText()
        object_box.blockSignals(True)
        object_box.clear()
        object_box.addItems(_eligible_objects())
        if current:
            idx = object_box.findText(current)
            if idx >= 0:
                object_box.setCurrentIndex(idx)
        object_box.blockSignals(False)
        backend_label.setText(f"当前编号后端：{_backend_cn(backend_status())}")

    def selected_color(combo, default):
        value = combo.currentData()
        return str(value) if value else default

    def save_colors():
        _prefs_set("cdr1_color", selected_color(c1, "yellow"))
        _prefs_set("cdr2_color", selected_color(c2, "magenta"))
        _prefs_set("cdr3_color", selected_color(c3, "cyan"))
        _prefs_set("show_sticks", sticks_cb.isChecked())

    def do_apply():
        save_colors()
        obj = object_box.currentText()
        if not obj:
            status.setPlainText("当前没有加载蛋白结构对象。")
            return
        try:
            result = _color_one(obj, quiet=True, show_sticks=sticks_cb.isChecked())
            status.setPlainText(format_result(result))
            _LAST_SIGNATURES[obj] = _object_signature(obj)
        except Exception as exc:
            status.setPlainText(f"识别失败：{exc}")

    def do_scan():
        save_colors()
        chunks = []
        for obj in _eligible_objects():
            try:
                result = _color_one(obj, quiet=True, show_sticks=sticks_cb.isChecked())
                chunks.append(format_result(result))
                _LAST_SIGNATURES[obj] = _object_signature(obj)
            except Exception as exc:
                chunks.append(f"{obj}：识别失败：{exc}")
        status.setPlainText("\n\n".join(chunks) if chunks else "当前没有加载蛋白结构对象。")
        refresh_objects()

    def do_clear():
        obj = object_box.currentText()
        if obj:
            delete_plugin_selections(obj)
            status.appendPlainText(f"已删除 {obj} 的插件 CDR selections；原子颜色不会自动恢复。")

    def toggle_auto(checked):
        cdr_auto(1 if checked else 0)

    apply_btn.clicked.connect(do_apply)
    scan_btn.clicked.connect(do_scan)
    clear_btn.clicked.connect(do_clear)
    auto_cb.toggled.connect(toggle_auto)
    sticks_cb.toggled.connect(lambda checked: _prefs_set("show_sticks", bool(checked)))

    timer = QtCore.QTimer(dialog)
    timer.setInterval(1500)
    timer.timeout.connect(refresh_objects)
    timer.start()
    dialog._refresh_timer = timer

    refresh_objects()
    _DIALOG = dialog
    dialog.show()


def __init_plugin__(app=None):
    """PyMOL 插件入口。"""
    from pymol import cmd
    from pymol.plugins import addmenuitemqt
    from .core import configure_hmmscan_path

    # macOS Finder 启动的 GUI 应用通常没有 Homebrew PATH。
    # 必须在自动 watcher 首次调用 ANARCI 之前完成注入。
    configure_hmmscan_path()

    cmd.extend("cdr_color", cdr_color)
    cmd.extend("cdr_info", cdr_info)
    cmd.extend("cdr_clear", cdr_clear)
    cmd.extend("cdr_auto", cdr_auto)
    cmd.extend("cdr_backend", cdr_backend)

    addmenuitemqt("抗体 CDR 自动识别与着色", run_plugin_gui)
    _start_watcher()
