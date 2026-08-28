"""抗体 CDR 自动识别与着色插件的核心逻辑。

默认采用 IMGT CDR 定义：
  CDR1: 27-38
  CDR2: 56-65
  CDR3: 105-117

编号后端均为可选依赖：优先 ANARCI，其次 ANARCII；如果二者都不可用，
则使用内置序列启发式算法。内置算法专门兼容 VHH/纳米抗体常见的长 CDR3
以及 CDR3 内额外 Cys，避免把额外 Cys 错当作 FR3 末端保守 Cys。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import os
import re
import shutil
import sys

IMGT_CDR_RANGES = {
    "CDR1": (27, 38),
    "CDR2": (56, 65),
    "CDR3": (105, 117),
}



# macOS 从 Finder 启动 PyMOL.app 时通常不会继承 Homebrew 的 PATH。
# ANARCI 内部通过命令名 ``hmmscan`` 调用 HMMER，因此即使终端里能找到
# /opt/homebrew/bin/hmmscan，GUI 版 PyMOL 仍可能报 FileNotFoundError。
_HMMSCAN_KNOWN_PATHS = (
    "/opt/homebrew/bin/hmmscan",   # Apple Silicon Homebrew
    "/usr/local/bin/hmmscan",     # Intel Homebrew
    "/opt/local/bin/hmmscan",     # MacPorts
    "/usr/bin/hmmscan",
)


def configure_hmmscan_path() -> Optional[str]:
    """寻找 hmmscan，并确保其目录位于当前 PyMOL 进程的 PATH。

    该函数只修改 PyMOL 当前进程的环境变量，不修改 shell 配置或系统文件。
    返回 hmmscan 的绝对路径；找不到则返回 None。
    """
    # 先尊重当前 PATH。
    found = shutil.which("hmmscan")
    if found:
        return os.path.realpath(found)

    candidates = []

    # PyMOL/当前 Python 环境自己的 bin 目录。
    for base in (
        os.path.dirname(sys.executable),
        os.path.join(sys.prefix, "bin"),
        os.path.join(getattr(sys, "base_prefix", sys.prefix), "bin"),
    ):
        if base:
            candidates.append(os.path.join(base, "hmmscan"))

    candidates.extend(_HMMSCAN_KNOWN_PATHS)

    # Conda/Mamba 环境有时通过这些变量暴露。
    for env_name in ("CONDA_PREFIX", "MAMBA_ROOT_PREFIX", "VIRTUAL_ENV"):
        prefix = os.environ.get(env_name)
        if prefix:
            candidates.append(os.path.join(prefix, "bin", "hmmscan"))

    seen = set()
    for path in candidates:
        path = os.path.expanduser(path)
        if path in seen:
            continue
        seen.add(path)
        if not (os.path.isfile(path) and os.access(path, os.X_OK)):
            continue

        bin_dir = os.path.dirname(path)
        current = os.environ.get("PATH", "")
        path_parts = [x for x in current.split(os.pathsep) if x]
        if bin_dir not in path_parts:
            os.environ["PATH"] = bin_dir + (os.pathsep + current if current else "")

        # 再通过 which 验证 ANARCI 用裸命令名启动时确实能解析到它。
        resolved = shutil.which("hmmscan")
        return os.path.realpath(resolved or path)

    return None


def hmmscan_status() -> Optional[str]:
    """返回当前 PyMOL 进程可用的 hmmscan 绝对路径。"""
    return configure_hmmscan_path()

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z",
}


@dataclass(frozen=True)
class ResidueRef:
    chain: str
    segi: str
    resi: str
    resn: str
    aa: str
    ca_index: int


@dataclass
class DomainResult:
    chain_key: Tuple[str, str]
    chain_type: str
    backend: str
    query_start: int
    query_end: int
    cdr_residues: Dict[str, List[ResidueRef]] = field(default_factory=dict)
    cdr_sequences: Dict[str, str] = field(default_factory=dict)


@dataclass
class ObjectResult:
    object_name: str
    backend: str
    domains: List[DomainResult]
    warnings: List[str] = field(default_factory=list)

    @property
    def has_antibody(self) -> bool:
        return bool(self.domains)

    @property
    def heavy_only(self) -> bool:
        """仅检测到重链型可变域；这类对象可能是 VH 或 VHH。"""
        if not self.domains:
            return False
        types = {d.chain_type.upper() for d in self.domains}
        return bool(types) and types <= {"H", "VHH", "VH"}


_ANARCII_MODEL = None


def backend_status() -> str:
    """返回当前真正可运行的最佳编号后端。

    ANARCI 只有在 Python 模块可导入且 hmmscan 可执行时才视为可用。
    """
    hmmscan = configure_hmmscan_path()
    if hmmscan:
        try:
            import anarci  # noqa: F401
            return "ANARCI"
        except Exception:
            pass
    try:
        import anarcii  # noqa: F401
        return "ANARCII"
    except Exception:
        return "heuristic"


def extract_chains(object_name: str, state: int = 1) -> Dict[Tuple[str, str], List[ResidueRef]]:
    """从 PyMOL 对象中按 CA 原子提取有序蛋白残基。"""
    from pymol import cmd

    model = cmd.get_model(f"%{object_name} and polymer.protein and name CA", state=state)
    chains: Dict[Tuple[str, str], List[ResidueRef]] = {}
    seen = set()

    for atom in model.atom:
        chain = getattr(atom, "chain", "") or ""
        segi = getattr(atom, "segi", "") or ""
        resi = str(getattr(atom, "resi", ""))
        key_res = (chain, segi, resi)
        if key_res in seen:
            continue
        seen.add(key_res)
        resn = str(getattr(atom, "resn", "UNK")).upper()
        aa = AA3_TO_1.get(resn, "X")
        ref = ResidueRef(
            chain=chain,
            segi=segi,
            resi=resi,
            resn=resn,
            aa=aa,
            ca_index=int(atom.index),
        )
        chains.setdefault((chain, segi), []).append(ref)

    return chains


def _loop_name(imgt_position: int) -> Optional[str]:
    for name, (lo, hi) in IMGT_CDR_RANGES.items():
        if lo <= imgt_position <= hi:
            return name
    return None


def _map_numbering_to_residues(
    residues: Sequence[ResidueRef],
    numbering: Sequence[Tuple[Tuple[int, str], str]],
    query_start: int,
    chain_type: str,
    backend: str,
) -> Optional[DomainResult]:
    if query_start is None or query_start < 0 or query_start >= len(residues):
        return None

    cdrs = {"CDR1": [], "CDR2": [], "CDR3": []}
    ptr = int(query_start)

    for item in numbering:
        if not item or len(item) != 2:
            continue
        pos_info, aa = item
        if aa == "-":
            continue
        if ptr >= len(residues):
            break
        try:
            pos = int(pos_info[0])
        except Exception:
            ptr += 1
            continue

        loop = _loop_name(pos)
        if loop:
            cdrs[loop].append(residues[ptr])
        ptr += 1

    if not any(cdrs.values()):
        return None

    end_idx = max(query_start, ptr - 1)
    return DomainResult(
        chain_key=(residues[0].chain, residues[0].segi),
        chain_type=(chain_type or "?").upper(),
        backend=backend,
        query_start=int(query_start),
        query_end=int(end_idx),
        cdr_residues=cdrs,
        cdr_sequences={k: "".join(r.aa for r in v) for k, v in cdrs.items()},
    )


def _number_with_anarci(sequence: str, residues: Sequence[ResidueRef]) -> List[DomainResult]:
    """使用经典 ANARCI 做 IMGT 编号；VHH 会按重链 H 处理。"""
    hmmscan = configure_hmmscan_path()
    if not hmmscan:
        raise FileNotFoundError(
            "未找到 hmmscan。已自动检查 PATH、/opt/homebrew/bin、"
            "/usr/local/bin、/opt/local/bin 及当前 Python/Conda 环境。"
        )

    from anarci import anarci

    # 不要求重/轻链成对出现；VHH 单链可直接作为 H 链编号。
    numbering, details, _ = anarci(
        [("pymol_chain", sequence)],
        scheme="imgt",
        output=False,
        allow={"H", "K", "L"},
    )
    if not numbering or numbering[0] is None:
        return []

    detail_list = details[0] if details and details[0] else []
    out: List[DomainResult] = []
    for i, domain in enumerate(numbering[0]):
        domain_numbering, start_idx, end_idx = domain
        detail = detail_list[i] if i < len(detail_list) else {}
        chain_type = detail.get("chain_type", "?") if isinstance(detail, dict) else "?"
        mapped = _map_numbering_to_residues(
            residues, domain_numbering, int(start_idx), chain_type, "ANARCI"
        )
        if mapped:
            mapped.query_end = int(end_idx)
            out.append(mapped)
    return out


def _get_anarcii_model():
    global _ANARCII_MODEL
    if _ANARCII_MODEL is None:
        from anarcii import Anarcii
        # 强制 CPU，避免 PyMOL 内部意外初始化 CUDA。
        _ANARCII_MODEL = Anarcii(
            seq_type="antibody", mode="accuracy", cpu=True, verbose=False
        )
    return _ANARCII_MODEL


def _number_with_anarcii(sequence: str, residues: Sequence[ResidueRef]) -> List[DomainResult]:
    model = _get_anarcii_model()
    results = model.number([sequence])
    if not isinstance(results, dict):
        return []

    out: List[DomainResult] = []
    for value in results.values():
        if not isinstance(value, dict):
            continue
        if value.get("error") or not value.get("numbering"):
            continue
        chain_type = str(value.get("chain_type", "?")).upper()
        if chain_type not in {"H", "K", "L"}:
            continue
        start_idx = value.get("query_start", 0)
        end_idx = value.get("query_end", len(sequence) - 1)
        mapped = _map_numbering_to_residues(
            residues,
            value["numbering"],
            int(start_idx or 0),
            chain_type,
            "ANARCII",
        )
        if mapped:
            mapped.query_end = int(end_idx if end_idx is not None else mapped.query_end)
            out.append(mapped)
    return out


def _find_fr4_motif(sequence: str):
    """寻找靠近可变域末端的 FR4 [F/W]GxG 样 motif。

    使用最后一个合理 motif，避免 VHH 较长 CDR3 内部偶然出现类似模式时被提前截断。
    """
    seq = sequence.upper()
    matches = [m for m in re.finditer(r"[FW]G[A-Z]G|[FW]GG", seq) if m.start() >= 75]
    if not matches:
        return None

    plausible = []
    for m in matches:
        tail_len = len(seq) - m.end()
        # FR4 后通常还会有少量残基；允许结构缺失或构建截短。
        if 0 <= tail_len <= 24:
            plausible.append(m)
    return (plausible or matches)[-1]


def _find_cdr3_anchor_cys(sequence: str, motif_start: int) -> Optional[int]:
    """寻找 FR3 末端的保守 Cys，而不是简单取 CDR3 前最后一个 Cys。

    VHH/纳米抗体 CDR3 内常存在额外 Cys，因此“取最后一个 Cys”会把 HCDR3
    错误截短。这里综合 Cys 所处区域、其前方 YY/Y 芳香族特征以及合理 CDR3
    长度评分，并在分数接近时偏向更靠前的 Cys。
    """
    seq = sequence.upper()
    candidates = []
    for i, aa in enumerate(seq[:motif_start]):
        if aa != "C" or i < 65:
            continue
        cdr3_len = motif_start - (i + 1)
        if not (1 <= cdr3_len <= 45):
            continue

        score = 0.0
        # FR3 末端保守 Cys 在未编号原始序列中通常位于约 85-110 区间。
        if 82 <= i <= 112:
            score += 3.0
        score -= abs(i - 100) * 0.04

        pre = seq[max(0, i - 6):i]
        if pre.endswith("YY"):
            score += 6.0
        elif pre.endswith("Y"):
            score += 4.0
        elif "Y" in pre[-3:]:
            score += 2.5
        if pre.endswith(("F", "W")):
            score += 1.0

        # 常见 VH/VHH HCDR3 长度优先，但仍给超长 VHH 留足空间。
        if 5 <= cdr3_len <= 30:
            score += 2.0
        if 8 <= cdr3_len <= 24:
            score += 1.0

        # tie-breaker：更靠前者更可能是保守 Cys，后面的 Cys 可能位于 CDR3 内。
        candidates.append((score, -i, i))

    if not candidates:
        return None
    candidates.sort(reverse=True)
    return candidates[0][2]


def _looks_like_variable_domain(sequence: str) -> bool:
    """保守的抗体可变域 motif gate，兼容 VHH 长 CDR3/额外 Cys。"""
    seq = sequence.upper()
    if len(seq) < 85:
        return False
    first_c = seq.find("C", 15, 40)
    if first_c < 0:
        return False

    aromatic_positions = [
        i for i in (
            seq.find("W", first_c + 4, first_c + 28),
            seq.find("Y", first_c + 4, first_c + 28),
        ) if i >= 0
    ]
    if not aromatic_positions:
        return False

    motif = _find_fr4_motif(seq)
    if motif is None:
        return False
    return _find_cdr3_anchor_cys(seq, motif.start()) is not None


def _guess_chain_type(sequence: str, chain_label: str) -> str:
    """根据序列优先判断重/轻链，链名只作为最后的弱提示。

    这样 VHH 即使在 PDB/mmCIF 中链 ID 恰好叫 L，也不会直接被误判成轻链。
    """
    seq = sequence.upper()

    # FR4 的首个芳香族通常是重链 W、轻链 F，是比任意 PDB 链名更可靠的提示。
    motif = _find_fr4_motif(seq)
    if motif is not None:
        marker = motif.group(0)[0]
        if marker == "W":
            return "H"
        if marker == "F":
            return "L"

    first_c = seq.find("C", 15, 40)
    window = seq[first_c:first_c + 60] if first_c >= 0 else seq[:75]
    if re.search(r"W[IVFLY][RK]Q|W[IVFLY]RQ|WFRQ|WVRQ|WIRQ|WVKQ", window):
        return "H"
    if "WYQ" in window or "WYLQ" in window:
        return "L"

    label = (chain_label or "").upper()
    if label in {"H", "VH", "VHH", "HC", "HEAVY"}:
        return "H"
    if label in {"L", "K", "VL", "LC", "LIGHT"}:
        return "L"
    return "?"


def _heuristic_cdrs(sequence: str, residues: Sequence[ResidueRef]) -> List[DomainResult]:
    """无依赖近似 CDR 划分，支持普通 Ig 和 VHH。

    该模式不是正式 IMGT 编号；正式分析建议安装 ANARCI/ANARCII。
    """
    seq = sequence.upper()
    if not _looks_like_variable_domain(seq):
        return []

    chain_type = _guess_chain_type(seq, residues[0].chain)
    if chain_type == "?":
        return []

    # 第一个保守 Cys（IMGT 23 附近）及其后的 FR2 芳香族 anchor。
    c1 = seq.find("C", 15, 40)
    aromatic = [
        i for i in (
            seq.find("W", c1 + 4, c1 + 28),
            seq.find("Y", c1 + 4, c1 + 28),
        ) if i >= 0
    ]
    if not aromatic:
        return []
    fr2_anchor = min(aromatic)

    # IMGT CDR1 从 27 开始；保守 Cys 通常为 23，因此跳过 Cys 后约 3 个 FR1 残基。
    cdr1_start = min(c1 + 4, fr2_anchor)
    cdr1_idx = list(range(cdr1_start, fr2_anchor))

    # IMGT CDR2 56-65。以 FR2 芳香族 anchor 近似定位。
    cdr2_start = fr2_anchor + 15
    cdr2_len = 10 if chain_type == "H" else 8
    cdr2_idx = list(range(cdr2_start, min(cdr2_start + cdr2_len, len(seq))))

    # CDR3：使用 FR4 motif + FR3 保守 Cys。关键点是允许 CDR3 内存在额外 Cys。
    motif_match = _find_fr4_motif(seq)
    if motif_match is None:
        return []
    c2 = _find_cdr3_anchor_cys(seq, motif_match.start())
    if c2 is None:
        return []
    cdr3_idx = list(range(c2 + 1, motif_match.start()))

    def refs(indices: Iterable[int]) -> List[ResidueRef]:
        return [residues[i] for i in indices if 0 <= i < len(residues)]

    cdrs = {
        "CDR1": refs(cdr1_idx),
        "CDR2": refs(cdr2_idx),
        "CDR3": refs(cdr3_idx),
    }
    if not cdrs["CDR1"] or not cdrs["CDR3"]:
        return []

    return [DomainResult(
        chain_key=(residues[0].chain, residues[0].segi),
        chain_type=chain_type,
        backend="heuristic",
        query_start=0,
        query_end=min(len(residues) - 1, motif_match.end() + 8),
        cdr_residues=cdrs,
        cdr_sequences={k: "".join(r.aa for r in v) for k, v in cdrs.items()},
    )]


def analyze_object(object_name: str, state: int = 1) -> ObjectResult:
    chains = extract_chains(object_name, state=state)
    domains: List[DomainResult] = []
    warnings: List[str] = []
    used_backends = set()

    for chain_key, residues in chains.items():
        if len(residues) < 50:
            continue
        seq = "".join(r.aa for r in residues)
        chain_domains: List[DomainResult] = []

        try:
            chain_domains = _number_with_anarci(seq, residues)
        except Exception as exc:
            warnings.append(
                f"链 {chain_key[0] or '<空链名>'}：ANARCI 不可用或编号失败：{exc}"
            )

        if not chain_domains:
            try:
                chain_domains = _number_with_anarcii(seq, residues)
            except Exception as exc:
                warnings.append(
                    f"链 {chain_key[0] or '<空链名>'}：ANARCII 不可用或编号失败：{exc}"
                )

        if not chain_domains:
            chain_domains = _heuristic_cdrs(seq, residues)

        domains.extend(chain_domains)
        used_backends.update(d.backend for d in chain_domains)

    backend = "+".join(sorted(used_backends)) if used_backends else backend_status()
    return ObjectResult(object_name=object_name, backend=backend, domains=domains, warnings=warnings)


def safe_name(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_]", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unnamed"


def _make_residue_selection(object_name: str, selection_name: str, ca_indices: Sequence[int]) -> int:
    """基于 CA atom index 创建整残基 selection。"""
    from pymol import cmd

    if not ca_indices:
        cmd.delete(selection_name)
        return 0

    tmp = cmd.get_unused_name("_ab_cdr_ca")
    try:
        cmd.select_list(tmp, object_name, list(map(int, ca_indices)), mode="index", quiet=1)
        count = cmd.select(selection_name, f"byres ?{tmp}", quiet=1)
    finally:
        cmd.delete(tmp)
    return int(count)


def _domain_role(chain_type: str) -> str:
    """将编号后端给出的链型归一化为 selection 使用的 H/L 标签。"""
    ct = (chain_type or "?").upper()
    if ct in {"H", "VH", "VHH"}:
        return "H"
    if ct in {"K", "L", "VL", "LIGHT"}:
        return "L"
    return ct if ct != "?" else "AB"


def apply_colors(
    result: ObjectResult,
    cdr1_color: str = "yellow",
    cdr2_color: str = "magenta",
    cdr3_color: str = "cyan",
    show_sticks: bool = False,
    selection_prefix: str = "abCDR",
) -> Dict[str, str]:
    """创建最少且不重复的 CDR selections，并完成着色。

    v1.2 起不再同时建立“对象汇总 CDR”与“分链 CDR”两套 selection：
    - 单个 H 型可变域（VHH/单独 VH）：仅 H_CDR1/2/3，共 3 个。
    - 常规 VH+VL：H_CDR1/2/3 + L_CDR1/2/3，共 6 个。
    - 同一对象存在多个同类型可变域时，自动加入原始 chain ID 以避免重名。
    """
    from pymol import cmd

    colors = {"CDR1": cdr1_color, "CDR2": cdr2_color, "CDR3": cdr3_color}
    obj_tag = safe_name(result.object_name)
    created: Dict[str, str] = {}

    # 升级/重复运行时先清掉本插件此前为该对象建立的 selection，
    # 防止 v1.0/v1.1 的“汇总 CDR”残留造成 VHH 看起来仍有 6 个 selection。
    delete_plugin_selections(result.object_name, prefix=selection_prefix)

    role_counts: Dict[str, int] = {}
    for domain in result.domains:
        role = _domain_role(domain.chain_type)
        role_counts[role] = role_counts.get(role, 0) + 1

    role_seen: Dict[str, int] = {}
    for domain in result.domains:
        role = _domain_role(domain.chain_type)
        role_seen[role] = role_seen.get(role, 0) + 1

        if role_counts.get(role, 0) <= 1:
            role_tag = role
        else:
            chain = safe_name(domain.chain_key[0] or domain.chain_key[1] or "blank")
            # 若同一 chain/segi 上确实有多个同型域，再加序号兜底。
            role_tag = f"{role}_{chain}"
            if role_seen[role] > 1 and any(
                key.startswith(f"{role_tag}_") for key in created
            ):
                role_tag = f"{role_tag}_{role_seen[role]}"

        for loop in ("CDR1", "CDR2", "CDR3"):
            refs = domain.cdr_residues.get(loop, [])
            if not refs:
                continue
            sel = safe_name(f"{selection_prefix}_{obj_tag}_{role_tag}_{loop}")
            _make_residue_selection(result.object_name, sel, [r.ca_index for r in refs])
            cmd.color(colors[loop], f"?{sel}")
            if show_sticks:
                cmd.show("sticks", f"?{sel}")
            created[f"{role_tag}_{loop}"] = sel

    return created


def delete_plugin_selections(object_name: Optional[str] = None, prefix: str = "abCDR") -> None:
    from pymol import cmd

    tag = safe_name(object_name) if object_name else None
    for name in cmd.get_names("selections"):
        if not name.startswith(prefix + "_"):
            continue
        if tag and not name.startswith(f"{prefix}_{tag}_"):
            continue
        cmd.delete(name)


def format_result(result: ObjectResult) -> str:
    if not result.domains:
        return f"{result.object_name}：未检测到抗体可变域（后端：{result.backend}）"

    format_note = "；仅检测到重链型可变域，兼容 VH/VHH 单域抗体" if result.heavy_only else ""
    lines = [
        f"{result.object_name}：检测到 {len(result.domains)} 个抗体可变域；后端：{result.backend}{format_note}"
    ]
    for i, d in enumerate(result.domains, 1):
        chain = d.chain_key[0] or "<空链名>"
        segi = d.chain_key[1]
        label = f"链 {chain}" + (f" / segi {segi}" if segi else "")
        type_text = "H（VH/VHH）" if d.chain_type.upper() == "H" else d.chain_type
        method_text = "内置启发式" if d.backend == "heuristic" else d.backend
        lines.append(f"  {i}. {label}：类型={type_text}，方法={method_text}")
        for loop in ("CDR1", "CDR2", "CDR3"):
            refs = d.cdr_residues.get(loop, [])
            if refs:
                lines.append(
                    f"     {loop}：{refs[0].resi}-{refs[-1].resi}  {d.cdr_sequences.get(loop, '')}"
                )
    if result.warnings:
        lines.append("  提示：")
        lines.extend(f"     - {w}" for w in result.warnings)
    return "\n".join(lines)
