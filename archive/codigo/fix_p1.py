with open("manuscript/SSEE_Paper1_Framework.tex", "r", encoding="utf-8") as f:
    content = f.read()

# Let's replace the equation with just text temporarily
content = content.replace("\\begin{equation}", "%\\begin{equation}")
content = content.replace("M_{\\rm dyn}(R) = \\Migimf(R)\\cdot\\KALz\\cdot\\left[1+f_\\nu\\right],", "%M_{\\rm dyn}(R) = \\Migimf(R)\\cdot\\KALz\\cdot\\left[1+f_\\nu\\right],")
content = content.replace("\\label{eq:mdyn}", "%\\label{eq:mdyn}")
content = content.replace("\\end{equation}", "%\\end{equation}")

with open("manuscript/SSEE_Paper1_Framework.tex", "w", encoding="utf-8") as f:
    f.write(content)
