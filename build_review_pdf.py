"""Build a polished submission-candidate PDF from the Paper-E1 manuscript."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path


PAPER = Path(__file__).resolve().parent
SOURCE = PAPER / "PAPER_E1_MANUSCRIPT_DRAFT_01_2026-08-09.md"
FIGURES = PAPER / "figures" / "output"
OUTDIR = PAPER / "output" / "pdf"
TMPDIR = PAPER / "tmp" / "pdfs"
TEX = TMPDIR / "paper_e1_submission_candidate.tex"
PDF = OUTDIR / "PAPER_E1_SUBMISSION_CANDIDATE_2026-08-09.pdf"

FIGURE_INSERTS = {
    "## 2. Observable-relative boundary-state framework": (
        "figure_1_method_architecture.pdf",
        "Extraction, predictive holdout, constrained boundary data, independent reconstruction, and geometric scoring.",
        "fig:architecture",
    ),
    "### 6.1 Predictive finite-window reconstruction": (
        "figure_2_production_convergence.pdf",
        "Common-surface depth and norm sensitivity. (a) The shallow depth-8 extraction fails while the deeper cases close under refinement; the depth-16 dashed and depth-24 solid curves coincide at displayed scale. (b) The conclusion is stable under identity, shape-emphasized, clock-emphasized, and component-maximum scores.",
        "fig:production",
    ),
    "### 6.2 Shape versus full marked geometry: the global modulus": (
        "figure_4_global_modulus_ablation.pdf",
        "Global-modulus ablation at (a) depth 16 and (b) depth 24. Shape denotes the maximum of the normalized component residuals; production denotes the full identity-weight score. Shape is unchanged, whereas marked-geometric closure fails when the anchor is removed.",
        "fig:ablation",
    ),
    "### 6.3 Exact Fuchsian and geometric benchmark": (
        "figure_3_metric_error_budget.pdf",
        "Exact Fuchsian and window/asymptotic reconstructions coincide at the discretization floor; the finite-slice metric route has a larger but second-order error.",
        "fig:error-budget",
    ),
    "## 7. Relation to Fuchsian and geometric singular-data methods": (
        "figure_5_modulus_recovery.pdf",
        "Modulus tests. (a) Known-answer estimator error follows the predicted asymptotic remainder. (b) Production estimates at depths 16 and 24 approach depth-consistent Richardson intercepts as inverse squared resolution tends to zero.",
        "fig:modulus",
    ),
    "### 6.4 Nonlinear unpolarized extension: localized insufficiency": (
        "figure_6_unpolarized_spike_extension.pdf",
        "Unpolarized extension. (a) The generic localized mismatch remains near 0.82 while the detected physical width is stable under refinement. (b) Ordinary four-function VTD data lose the exact spike, whereas the enriched centered-spike state restores the tested exact family.",
        "fig:unpolarized-spike",
    ),
}


def escape_text(text: str) -> str:
    protected = []

    def protect(pattern, replacement):
        nonlocal text
        def repl(match):
            protected.append(replacement(match))
            return f"@@PROTECTED{len(protected)-1}@@"
        text = re.sub(pattern, repl, text)

    protect(r"\\\((.+?)\\\)", lambda m: f"${m.group(1)}$")
    protect(r"`([^`]+)`", lambda m: r"\texttt{" + latex_escape(m.group(1)) + "}")
    protect(r"\[([^\]]+)\]\(([^)]+)\)",
            lambda m: r"\href{" + m.group(2) + "}{" + latex_escape(m.group(1)) + "}")
    text = latex_escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", text)
    text = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"\\emph{\1}", text)
    for i, value in enumerate(protected):
        text = text.replace(f"@@PROTECTED{i}@@", value)
    return text


def latex_escape(text: str) -> str:
    mapping = {
        "&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_",
        "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(mapping.get(char, char) for char in text)


def table_to_latex(lines: list[str]) -> str:
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in lines]
    if len(rows) > 1 and all(re.fullmatch(r":?-+:?", c.replace(" ", "")) for c in rows[1]):
        rows.pop(1)
    cols = len(rows[0])
    spec = "@{}" + "l" + "Y" * (cols - 1) + "@{}"
    out = [r"\begin{table}[H]", r"\centering", r"\small",
           rf"\begin{{tabularx}}{{\textwidth}}{{{spec}}}", r"\toprule"]
    for idx, row in enumerate(rows):
        cells = [escape_text(cell) for cell in row]
        if idx == 0:
            cells = [r"\textbf{" + c + "}" for c in cells]
        out.append(" & ".join(cells) + r" \\")
        if idx == 0:
            out.append(r"\midrule")
    out.extend([r"\bottomrule", r"\end{tabularx}", r"\end{table}"])
    return "\n".join(out)


def figure_block(filename: str, caption: str, label: str) -> str:
    path = (FIGURES / filename).as_posix()
    return "\n".join([
        r"\begin{figure}[H]", r"\centering",
        rf"\includegraphics[width=0.97\textwidth]{{{path}}}",
        rf"\caption{{{latex_escape(caption)}}}", rf"\label{{{label}}}",
        r"\end{figure}",
    ])


def convert(markdown: str) -> str:
    lines = markdown.splitlines()
    out: list[str] = []
    paragraph: list[str] = []
    in_math = False
    in_code = False
    list_kind = None

    def flush_paragraph():
        if paragraph:
            out.append(escape_text(" ".join(s.strip() for s in paragraph)))
            out.append("")
            paragraph.clear()

    def close_list():
        nonlocal list_kind
        if list_kind:
            out.append(rf"\end{{{list_kind}}}")
            list_kind = None

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if stripped.startswith(("**Manuscript status:**", "**Provisional target:**",
                                "**Target:**", "**Article class:**", "**Author:**",
                                "**Authors:**", "**Affiliation:**", "**Email:**",
                                "**ORCID:**")):
            i += 1
            continue

        if stripped.startswith("```"):
            flush_paragraph(); close_list()
            if not in_code:
                out.append(r"\begin{Verbatim}[fontsize=\small]")
            else:
                out.append(r"\end{Verbatim}")
            in_code = not in_code
            i += 1
            continue
        if in_code:
            out.append(line)
            i += 1
            continue

        if stripped == r"\[":
            flush_paragraph(); close_list(); in_math = True; out.append(r"\[")
            i += 1; continue
        if stripped == r"\]":
            out.append(r"\]"); in_math = False; i += 1; continue
        if in_math:
            out.append(line)
            i += 1; continue

        if stripped.startswith("|") and i + 1 < len(lines) and lines[i + 1].strip().startswith("|"):
            flush_paragraph(); close_list()
            table_lines = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                table_lines.append(lines[i]); i += 1
            out.append(table_to_latex(table_lines)); continue

        if stripped in FIGURE_INSERTS:
            flush_paragraph(); close_list()
            out.append(figure_block(*FIGURE_INSERTS[stripped]))

        if stripped == "## Draft completion ledger":
            flush_paragraph(); close_list()
            break
        if stripped.startswith("# "):
            flush_paragraph(); close_list()
            out.append(r"\title{" + escape_text(stripped[2:]) + "}")
        elif stripped.startswith("## "):
            flush_paragraph(); close_list()
            title = re.sub(r"^\d+\.\s*", "", stripped[3:])
            appendix_match = re.match(r"^Appendix ([A-Z])\.\s*(.*)$", title)
            if appendix_match:
                if appendix_match.group(1) == "A":
                    out.append(r"\appendix")
                title = appendix_match.group(2)
            if title in ("Abstract", "Acknowledgments", "Data and code availability", "References"):
                out.append(r"\section*{" + escape_text(title) + "}")
            else:
                out.append(r"\section{" + escape_text(title) + "}")
        elif stripped.startswith("### "):
            flush_paragraph(); close_list()
            title = re.sub(r"^\d+(?:\.\d+)*\s*", "", stripped[4:])
            title = re.sub(r"^[A-Z]\.\d+\s*", "", title)
            out.append(r"\subsection{" + escape_text(title) + "}")
        elif re.match(r"^- \[[ xX]\]", stripped):
            flush_paragraph()
            if list_kind != "itemize": close_list(); out.append(r"\begin{itemize}"); list_kind = "itemize"
            checked = stripped[3].lower() == "x"
            out.append(r"\item[" + (r"$\checkmark$" if checked else r"$\square$") + "] " + escape_text(stripped[6:]))
        elif re.match(r"^- ", stripped):
            flush_paragraph()
            if list_kind != "itemize": close_list(); out.append(r"\begin{itemize}"); list_kind = "itemize"
            out.append(r"\item " + escape_text(stripped[2:]))
        elif re.match(r"^\d+\. ", stripped):
            flush_paragraph()
            if list_kind != "enumerate": close_list(); out.append(r"\begin{enumerate}"); list_kind = "enumerate"
            out.append(r"\item " + escape_text(re.sub(r"^\d+\. ", "", stripped)))
        elif stripped == "":
            flush_paragraph(); close_list()
        elif stripped.startswith("> "):
            flush_paragraph(); close_list()
            out.append(r"\begin{quote}" + escape_text(stripped[2:]) + r"\end{quote}")
        else:
            paragraph.append(stripped)
        i += 1

    flush_paragraph(); close_list()
    return "\n".join(out)


PREAMBLE = r"""
\documentclass[11pt]{article}
\usepackage[margin=0.85in]{geometry}
\usepackage{fontspec}
\setmainfont{STIXGeneral}
\setsansfont{Arial}
\setmonofont{Courier New}
\usepackage{amsmath,amssymb,mathtools}
\usepackage{unicode-math}
\setmathfont{STIX Two Math}
\usepackage{graphicx,float,booktabs,tabularx,array,longtable}
\usepackage{xcolor,hyperref,fancyhdr,microtype,enumitem,fancyvrb}
\usepackage{pifont}
\newcolumntype{Y}{>{\raggedright\arraybackslash}X}
\hypersetup{colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black,citecolor=blue!50!black}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{Boundary-State Sufficiency in Gowdy Cosmology}
\fancyhead[R]{Ronald Bibb}
\fancyfoot[C]{\thepage}
\setlength{\headheight}{14pt}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\setlength{\emergencystretch}{2em}
\setlist{nosep,leftmargin=1.6em}
\author{Ronald Bibb\\
\normalsize Independent Researcher, Lilburn, GA 30047, USA\\
\normalsize \href{mailto:ronbibb@gmail.com}{ronbibb@gmail.com}\\
\normalsize ORCID: \href{https://orcid.org/0009-0004-1153-2464}{0009-0004-1153-2464}}
\date{}
\begin{document}
\maketitle
\vspace{8pt}
"""

POSTAMBLE = r"\end{document}"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)
    body = convert(SOURCE.read_text())
    # Markdown title is converted to \title and must precede \maketitle.
    title_match = re.search(r"\\title\{.*?\}\n", body)
    title = title_match.group(0) if title_match else r"\title{Paper E1}\n"
    body = body.replace(title, "", 1)
    tex = PREAMBLE.replace(r"\begin{document}", title + r"\begin{document}", 1) + body + "\n" + POSTAMBLE
    TEX.write_text(tex)
    for _ in range(2):
        subprocess.run(
            ["/Library/TeX/texbin/xelatex", "-interaction=nonstopmode", "-halt-on-error", TEX.name],
            cwd=TMPDIR, check=True,
        )
    built = TMPDIR / TEX.with_suffix(".pdf").name
    PDF.write_bytes(built.read_bytes())
    print(PDF)


if __name__ == "__main__":
    main()
