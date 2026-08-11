"""Build an Overleaf-ready IOP/CQG LaTeX source from the E1 manuscript."""

from __future__ import annotations

import re
from pathlib import Path

from build_review_pdf import PAPER, SOURCE, convert


OUTPUT = PAPER / "CQG_OVERLEAF_SUBMISSION_2026-08-10.tex"


def extract_section(body: str, heading: str, next_heading: str) -> tuple[str, str]:
    pattern = (
        rf"\\section\*\{{{re.escape(heading)}\}}\n"
        rf"(.*?)"
        rf"(?=\\(?:section\*?|appendix)\{{{re.escape(next_heading)}\}}|"
        rf"\\section\{{{re.escape(next_heading)}\}})"
    )
    match = re.search(pattern, body, flags=re.S)
    if not match:
        raise ValueError(f"Could not extract {heading!r}")
    return match.group(1).strip(), body[:match.start()] + body[match.end():]


def main() -> None:
    body = convert(SOURCE.read_text())

    title_match = re.search(r"\\title\{(.*?)\}\n", body)
    if not title_match:
        raise ValueError("Manuscript title not found")
    title = title_match.group(1)
    body = body[:title_match.start()] + body[title_match.end():]

    abstract_match = re.search(
        r"\\section\*\{Abstract\}\n(.*?)(?=\\section\{Introduction\})",
        body, flags=re.S,
    )
    if not abstract_match:
        raise ValueError("Abstract not found")
    abstract = abstract_match.group(1).strip()
    body = body[:abstract_match.start()] + body[abstract_match.end():]

    ack_match = re.search(
        r"\\section\*\{Acknowledgments\}\n(.*?)"
        r"(?=\\section\*\{Data and code availability\})",
        body, flags=re.S,
    )
    if not ack_match:
        raise ValueError("Acknowledgments not found")
    acknowledgments = ack_match.group(1).strip()
    body = body[:ack_match.start()] + body[ack_match.end():]

    data_match = re.search(
        r"\\section\*\{Data and code availability\}\n(.*?)"
        r"(?=\\section\*\{References\})",
        body, flags=re.S,
    )
    if not data_match:
        raise ValueError("Data statement not found")
    data_statement = data_match.group(1).strip()
    body = body[:data_match.start()] + body[data_match.end():]

    # The review builder uses absolute figure paths for local rendering.
    # Overleaf receives the repository tree, so journal source must be relative.
    body = body.replace(PAPER.as_posix() + "/", "")

    references_marker = r"\section*{References}"
    body = body.replace(
        references_marker,
        "\n".join([
            rf"\ack{{{acknowledgments}}}",
            r"\funding{This research received no institutional or grant funding.}",
            rf"\data{{{data_statement}}}",
            "",
            references_marker,
        ]),
        1,
    )

    preamble = rf"""% Overleaf main document for submission to Classical and Quantum Gravity.
% Uses IOP Publishing's official iopjournal class (2025 template package).
\documentclass{{iopjournal}}

\usepackage{{amsmath,amssymb,mathtools,mathrsfs}}
\usepackage{{booktabs,tabularx,array,longtable,float}}
\usepackage{{enumitem,fancyvrb,microtype,pifont}}
\newcolumntype{{Y}}{{>{{\raggedright\arraybackslash}}X}}
\setlength{{\emergencystretch}}{{2em}}
\setlist{{nosep,leftmargin=1.6em}}
\hypersetup{{pdftitle={{{title}}},pdfauthor={{Ronald Bibb}}}}

\begin{{document}}

\articletype{{Paper}}

\title{{{title}}}

\author{{Ronald Bibb\orcid{{0009-0004-1153-2464}}}}
\affil{{Independent Researcher, Lilburn, GA 30047, USA}}
\email{{ronbibb@gmail.com}}
\par\noindent\keywords{{cosmological singularities, Gowdy spacetimes, Fuchsian asymptotics, numerical relativity, spike solutions}}\par

\begin{{abstract}}
{abstract}
\end{{abstract}}

"""
    OUTPUT.write_text(preamble + body.strip() + "\n\n\\end{document}\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
