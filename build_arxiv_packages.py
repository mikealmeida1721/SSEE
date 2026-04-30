import os
import shutil
import re
import tarfile

papers = [
    {"id": 1, "name": "SSEE_Paper1_Framework", "files": ["SSEE_Paper1_Framework.tex", "SSEE_EFT_section.tex"]},
    {"id": 2, "name": "SSEE_Paper2_MCMC", "files": ["SSEE_Paper2_MCMC.tex", "SSEE_Paper2_MCMC.bbl"]},
    {"id": 3, "name": "SSEE_Paper3_CMB", "files": ["SSEE_Paper3_CMB.tex", "SSEE_Paper3_CMB.bbl"]},
    {"id": 4, "name": "SSEE_Paper4_ToE", "files": ["SSEE_Paper4_ToE.tex", "SSEE_Paper4_ToE.bbl"]}
]

os.makedirs("submission_packages", exist_ok=True)

for p in papers:
    pid = p["id"]
    pname = p["name"]
    pkg_dir = f"pkg_{pname}"
    os.makedirs(pkg_dir, exist_ok=True)
    os.makedirs(f"{pkg_dir}/figures", exist_ok=True)
    
    # Read main tex file
    main_tex_path = f"manuscript/{pname}.tex"
    if not os.path.exists(main_tex_path):
        continue
        
    with open(main_tex_path, "r") as f:
        content = f.read()
        
    # Find all includegraphics
    figs = re.findall(r'\\includegraphics(?:\[.*?\])?\{([^}]+)\}', content)
    
    # Process figures
    for fig in figs:
        # Resolve fig path
        fig_clean = fig.replace("../results/figures/", "").replace("figures/", "")
        src_fig = f"results/figures/{fig_clean}"
        dst_fig = f"{pkg_dir}/figures/{fig_clean}"
        if os.path.exists(src_fig):
            shutil.copy(src_fig, dst_fig)
            
    # Rewrite paths in tex content
    content = content.replace("../results/figures/", "figures/")
    content = re.sub(r'\\graphicspath\{\{.*?\}\}', r'\\graphicspath{{figures/}}', content)
    
    # Copy and rewrite tex files
    for fname in p["files"]:
        if fname == f"{pname}.tex":
            with open(f"{pkg_dir}/{fname}", "w") as f:
                f.write(content)
        else:
            src = f"manuscript/{fname}"
            if os.path.exists(src):
                shutil.copy(src, f"{pkg_dir}/{fname}")
                
    # Create tar.gz
    tar_name = f"submission_packages/{pname}_arXiv.tar.gz"
    with tarfile.open(tar_name, "w:gz") as tar:
        for item in os.listdir(pkg_dir):
            tar.add(f"{pkg_dir}/{item}", arcname=item)
            
    shutil.rmtree(pkg_dir)

print("Submission packages created successfully.")
