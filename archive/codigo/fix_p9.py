with open("manuscript/SSEE_Paper9_HubbleTension.tex", "r", encoding="utf-8") as f:
    content = f.read()

# Add \usepackage{tabularx}
content = content.replace("\\usepackage{graphicx}", "\\usepackage{graphicx}\n\\usepackage{tabularx}")

# Add \sloppy
if "\\sloppy" not in content:
    content = content.replace("\\begin{document}", "\\sloppy\n\\begin{document}")

# Fix specific table
import re
pattern = r"\\begin\{tabular\}\{lp\{5cm\}l\}(.*?)\\end\{tabular\}"
replacement = r"\\begin{tabularx}{\\textwidth}{l p{4.5cm} X}\1\\end{tabularx}"
content = re.sub(pattern, replacement, content, flags=re.DOTALL)

with open("manuscript/SSEE_Paper9_HubbleTension.tex", "w", encoding="utf-8") as f:
    f.write(content)
