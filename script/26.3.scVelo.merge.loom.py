# merge loompy  
import loompy
import glob
from pathlib import Path
import os


dir_main = Path("/share/home/shh/HCC_snRNA-seq")
# dir_in = dir_main.joinpath("26.velocyto/dataset")
dir_dataset = dir_main.joinpath("26.scVelo/dataset")
dir_result = dir_main.joinpath("26.scVelo/result")
os.chdir(str(dir_dataset))


loom_files = glob.glob(str(dir_dataset.joinpath('*.loom')))
loom_files = [file for file in loom_files if file.find("merge") == -1]
loompy.combine(loom_files, "merge.loom", key="Accession")
