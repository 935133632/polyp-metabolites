"""Assemble the publication package in publication_package/ (manuscript .docx, tables, figures, supplementary information, references, code snapshot, checklist)."""
import sys, re, csv, json, shutil, zipfile, importlib.util
sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
from pathlib import Path
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_COLOR_INDEX, WD_BREAK
from docx.shared import Cm, Pt
import pandas as pd
from docx_helpers import *

R = Path(r"D:\research\Xingang Li\polyp_metabolites"); PKG = R / "publication_package"
SP = Path(r"C:\Users\93513\AppData\Local\Temp\claude\d--research-Xingang-Li-polyp-metabolites\50704e06-a4d0-4ac0-83f3-1c38f08ea861\scratchpad")
spec = importlib.util.spec_from_file_location("main_tables", R / "code" / "19_main_tables.py"); mt = importlib.util.module_from_spec(spec); spec.loader.exec_module(mt)
LAB = mt.LAB
for d in ("01_Manuscript", "02_Figures/Main", "03_Tables", "04_Supplementary/Supplementary_Figures", "05_References", "06_Code_and_results"): (PKG / d).mkdir(parents=True, exist_ok=True)
md = (R / "manuscript" / "polyp_metabolites_draft.md").read_text(encoding="utf8"); legends = (R / "manuscript" / "figure_legends.md").read_text(encoding="utf8")
main_leg, supp_leg = legends.split("## Supplementary figures")
leg_items = lambda t: [x.strip() for x in re.split(r"\n\s*\n", t) if x.strip().startswith("**Figure")]

DECL = [("Ethics approval and consent to participate", "UK Biobank received ethical approval from the North West Multi-centre Research Ethics Committee, and all participants gave written informed consent. This research was conducted under UK Biobank application number [to be completed]."),
        ("Data availability", "UK Biobank data are available to approved researchers through the UK Biobank access procedure (www.ukbiobank.ac.uk). Aggregate results are provided in the Supplementary Tables."),
        ("Code availability", "Analysis code is available at https://github.com/935133632/polyp-metabolites."),
        ("Funding", "[to be completed]"), ("Competing interests", "[to be completed]"), ("Author contributions", "[to be completed]"),
        ("Acknowledgements", "We thank the UK Biobank participants and the UK Biobank team.")]

def para(doc, text, size=None, hanging=False, align=None, keep=False):
    p = doc.add_paragraph()
    for tok in re.split(r"(\[to be completed\])", text):
        if tok == "[to be completed]": r = p.add_run(tok); r.font.highlight_color = WD_COLOR_INDEX.YELLOW
        elif tok: add_runs(p, tok, size=size)
    if hanging: p.paragraph_format.left_indent = Cm(0.8); p.paragraph_format.first_line_indent = Cm(-0.8); p.paragraph_format.line_spacing = 1.15
    if align: p.alignment = align
    return p

# ----------------------------------------------------------------------------- 1. main manuscript
doc = new_document(size=11, line_spacing=1.5)
title = re.match(r"# (.+)", md).group(1)
t = doc.add_paragraph(); t.alignment = WD_ALIGN_PARAGRAPH.CENTER; r = t.add_run(title); r.bold = True; r.font.size = Pt(15)
para(doc, "[Author names, affiliations and corresponding author: to be completed]", align=WD_ALIGN_PARAGRAPH.CENTER)
body = md[md.index("## Abstract"):]; main_text, refs_text = body.split("## References")
for block in re.split(r"\n\s*\n", main_text.strip()):
    lines = block.strip().split("\n"); first = lines[0]          # a heading line may be followed directly by its paragraph
    if first.startswith("## "): doc.add_heading(first[3:], level=1); lines = lines[1:]
    elif first.startswith("### "): doc.add_heading(first[4:], level=2); lines = lines[1:]
    if lines and "".join(lines).strip(): para(doc, " ".join(l.strip() for l in lines))
doc.add_heading("Declarations", level=1)
for h, txt in DECL:
    p = doc.add_paragraph(); r = p.add_run(h + ". "); r.bold = True
    for tok in re.split(r"(\[to be completed\])", txt):
        if tok == "[to be completed]": rr = p.add_run(tok); rr.font.highlight_color = WD_COLOR_INDEX.YELLOW
        elif tok: p.add_run(tok)
doc.add_heading("References", level=1)
n_refs = 0
for line in refs_text.strip().splitlines():
    if line.strip(): para(doc, line.strip(), size=10, hanging=True); n_refs += 1
# tables (landscape)
sec = doc.add_section(WD_SECTION.NEW_PAGE); landscape(sec)
VARS, cols, rows = mt.table1()
head1 = ["Characteristic"] + [LAB["K63"], "", "", "", LAB["K31"], "", "", ""]
head2 = [""] + sum([[f"All participants\n(n = {cols[k]['All']:,})", f"Incident polyp\n(n = {cols[k]['Case']:,})", f"No polyp\n(n = {cols[k]['None']:,})", "P"] for k in LAB], [])
body1 = [[lab] + rows[lab]["K63"] + rows[lab]["K31"] for lab, *_ in VARS]
indent = [i for i, v in enumerate(VARS) if v[2].startswith(("cntNever", "cntPrevious", "cntCurrent"))]; heads = [i for i, v in enumerate(VARS) if v[2] == "head"]
caption(doc, "**Table 1.** Baseline characteristics of UK Biobank participants by incident colon polyp and incident stomach/duodenal polyp status", size=10)
three_line_table(doc, [head1, head2], body1, col_widths_cm=[6.2] + [2.2, 2.2, 2.2, 1.3] * 2, font_size=7.5, merges=[(0, 1, 4), (0, 5, 8)], indent_rows=indent, group_rows=heads)
caption(doc, "Values are mean (standard deviation), median (interquartile range) or n (%); percentages exclude missing values. Participants with the outcome before baseline were excluded from the corresponding cohort. P values compare participants with and without incident polyp (Welch t test, Mann–Whitney test or chi-squared test). Prior hospital diagnoses are those recorded before baseline. Metabolite concentrations are NMR measurements.", size=8)
doc.add_page_break()
body2, bold2, groups2, ns = mt.table2()
h1 = ["Metabolite", LAB["K63"], "", "", "", LAB["K31"], "", "", ""]; h2 = [""] + ["Model 1", "Model 2", "Model 3", "Model 4"] * 2
caption(doc, "**Table 2.** Associations of key metabolites with incident colon polyp and stomach/duodenal polyp across adjustment models (hazard ratio per 1 SD, 95% CI)", size=10)
three_line_table(doc, [h1, h2], body2, col_widths_cm=[6.0] + [2.5] * 8, font_size=7.5, merges=[(0, 1, 4), (0, 5, 8)], bold_cells=bold2, indent_rows=[i for i in range(len(body2)) if i not in groups2], group_rows=groups2)
caption(doc, f"Same complete-case sample for all models (colon polyp: {ns['K63'][0]:,} participants, {ns['K63'][1]:,} events; stomach/duodenal polyp: {ns['K31'][0]:,}, {ns['K31'][1]:,}). Model 1: age, sex, Townsend deprivation index, body mass index, systolic and diastolic blood pressure, smoking and alcohol drinking. Model 2: Model 1 plus waist circumference, body fat percentage, HbA1c and C-reactive protein. Model 3: Model 1 plus prior hospital diagnoses (log number of ICD-10 codes; disorders of lipoprotein metabolism, ischaemic heart disease, hypertension, diabetes, reflux disease, gastritis or ulcer). Model 4: Models 2 and 3 combined. Metabolites were rank-normalised and standardised. Bold: P < 0.05 after Bonferroni correction across 251 metabolites. Results for all 251 metabolites are in Table S1 and Figure 3.", size=8)
# figure legends (portrait)
sec = doc.add_section(WD_SECTION.NEW_PAGE); sec.orientation = WD_ORIENT.PORTRAIT; sec.page_width, sec.page_height = sec.page_height, sec.page_width; sec.left_margin = sec.right_margin = Cm(2.5)
doc.add_heading("Figure legends", level=1)
for item in leg_items(main_leg): para(doc, item, size=10)
doc.save(PKG / "01_Manuscript" / "Manuscript.docx")
# tables as separate files
for f in ("Table1_baseline_characteristics.docx", "Table2_key_metabolites.docx", "Main_Tables.xlsx"): shutil.copy(R / "tables" / f, PKG / "03_Tables" / f)

# ----------------------------------------------------------------------------- 2. figures
FIG = R / "figures"; SUP = FIG / "supplement"
main_files = sorted(FIG.glob("Fig*.png"))
for f in main_files:
    for ext in (".png", ".pdf"):
        if f.with_suffix(ext).exists(): shutil.copy(f.with_suffix(ext), PKG / "02_Figures/Main" / f.with_suffix(ext).name)
for f in SUP.glob("FigS*"): shutil.copy(f, PKG / "04_Supplementary/Supplementary_Figures" / f.name)

# ----------------------------------------------------------------------------- 3. supplementary information
sd = new_document(size=11, line_spacing=1.3)
t = sd.add_paragraph(); t.alignment = WD_ALIGN_PARAGRAPH.CENTER; r = t.add_run("Supplementary Information"); r.bold = True; r.font.size = Pt(15)
para(sd, f"*{title}*", align=WD_ALIGN_PARAGRAPH.CENTER)
sd.add_heading("Supplementary tables", level=1)
para(sd, "Supplementary tables are provided in the Excel workbook **Supplementary_Tables.xlsx** (one sheet per table, with an index sheet).")
idx = pd.read_excel(R / "tables" / "supplement" / "Supplementary_Tables.xlsx", sheet_name="Index")
three_line_table(sd, [["Sheet", "Content"]], [[a, b] for a, b in zip(idx.sheet, idx.description)], col_widths_cm=[5.0, 11.0], font_size=9)
sd.add_page_break(); sd.add_heading("Supplementary figures", level=1)
legs = leg_items(supp_leg)
FILES = {"S1": ["FigS1A*.png", "FigS1B*.png"], "S2": ["FigS2A*.png", "FigS2B*.png"], "S3": ["FigS3*.png"], "S4": ["FigS4*.png"], "S5": ["FigS5*.png"], "S6": ["FigS6*.png"], "S7": ["FigS7*.png"], "S8": ["FigS8_*K63.png", "FigS8_*K31.png"], "S9": ["FigS9*.png"]}
for leg in legs:
    key = re.match(r"\*\*Figure (S\d+)\.", leg).group(1)
    for pat in FILES[key]:
        for img in sorted(SUP.glob(pat)): sd.add_picture(str(img), width=Cm(15 if key not in ("S4", "S5", "S9") else 15.5)); sd.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    para(sd, leg, size=10); sd.add_page_break()
sd.save(PKG / "04_Supplementary" / "Supplementary_Information.docx")
shutil.copy(R / "tables" / "supplement" / "Supplementary_Tables.xlsx", PKG / "04_Supplementary" / "Supplementary_Tables.xlsx")

# ----------------------------------------------------------------------------- 4. references (RIS + csv, in citation order)
refs = json.loads((SP / "refs_final.json").read_text(encoding="utf8")); order = list(csv.DictReader(open(R / "manuscript" / "references_verified.csv", encoding="utf8")))
doi2key = {v["doi"].lower(): k for k, v in refs.items()}; lines = []
for row in order:
    r_ = refs[doi2key[row["doi"].lower()]]; ty = "BOOK" if r_["key"] == "harrell2015" else "CPAPER" if r_["key"] == "xgboost2016" else "JOUR"
    lines += [f"TY  - {ty}"] + [f"AU  - {a['family'].strip(' ,').title() if a['family'].isupper() else a['family'].strip(' ,')}" + (f", {a['given']}" if a["given"] else "") for a in r_["authors"]] + [f"TI  - {r_['title']}", f"{'T2' if ty != 'JOUR' else 'JO'}  - {r_['journal'] or r_.get('publisher') or ''}"]
    if r_.get("volume"): lines.append(f"VL  - {r_['volume']}")
    if r_.get("issue"): lines.append(f"IS  - {r_['issue']}")
    pg = (r_.get("page") or "").split("-")
    if pg and pg[0]: lines.append(f"SP  - {pg[0]}")
    if len(pg) > 1: lines.append(f"EP  - {pg[1]}")
    lines += [f"PY  - {r_['year']}", f"DO  - {r_['doi']}", f"UR  - https://doi.org/{r_['doi']}", "ER  - ", ""]
(PKG / "05_References" / "references.ris").write_text("\n".join(lines), encoding="utf8"); shutil.copy(R / "manuscript" / "references_verified.csv", PKG / "05_References" / "references.csv")

# ----------------------------------------------------------------------------- 5. code and aggregate results snapshot
z = PKG / "06_Code_and_results" / "code_and_aggregate_results.zip"
with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
    for f in sorted((R / "code").glob("*")):
        if f.is_file() and f.suffix in (".py", ".ps1"): zf.write(f, f"code/{f.name}")
    for f in sorted((R / "results").rglob("*.csv")):
        if "bootstrap_chunks" not in f.parts: zf.write(f, f"results/{f.relative_to(R / 'results').as_posix()}")
    zf.write(R / "ANALYSIS_LOG.md", "ANALYSIS_LOG.md"); zf.write(R / "README.md", "README.md")
print("package built:", sum(1 for _ in PKG.rglob("*") if _.is_file()), "files |", n_refs, "references |", len(main_files), "main figure images | zip", round(z.stat().st_size / 1e6, 1), "MB")
