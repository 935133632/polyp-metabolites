"""Rule-based class annotation of the 251 Nightingale NMR metabolite abbreviations used in this dataset."""
import re
AA = {"Alanine", "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Valine", "Phenylalanine", "Tyrosine", "Concen_BCAA"}
FA = {"FattyAcids", "Unsaturation", "Omega3", "Omega6", "Polyunsaturated", "Monounsaturated", "Saturated", "Linoleic", "Docosahexaenoic",
      "O3pt", "O6pt", "Polypt", "Monopt", "Satupt", "Linopt", "Docopt", "PolyRatio", "O6Ratio"}
GLY = {"Glucose", "Lactate", "Pyruvate", "Citrate", "3Hydroxybutyrate", "Acetate", "Acetoacetate", "Acetone", "GluLac"}
APO = {"ApoB", "ApoA1", "ApoRatio"}
MISC = {"Creatinine": "Fluid balance", "Albumin": "Fluid balance", "GlycoAcet": "Inflammation (GlycA)"}


def metab_class(n: str) -> str:
    if n in MISC: return MISC[n]
    if n in AA: return "Amino acids"
    if n in FA: return "Fatty acids"
    if n in GLY: return "Glycolysis / ketone bodies"
    if n in APO: return "Apolipoproteins"
    if n.startswith("Diame_"): return "Lipoprotein particle size"
    if n.startswith("Concen_"): return "Lipoprotein particle concentration"
    if re.search(r"_pt$", n): return "Lipid composition (% within subclass)"
    if re.match(r"^(Lip_|Lip__|Pho_|Cho_|CE_|FC_|TG_|TC_)\w+", n) and re.search(r"(VLDL|LDL|HDL|IDL|Chy)", n): return "Lipid content of subclass"
    return "Lipoprotein / lipid totals"


def lipo_family(n: str) -> str:
    for f in ("VLDL", "IDL", "LDL", "HDL"):
        if f in n: return f
    return "Other"
