import os

rootdir = "/home/cannons/interpretation_graphs/consurf_scores"

for subdir, dir, files in os.walk(rootdir):
    for f in files:
        if f == "consurf.grades" or f == "keep_consurf.py":
            print(f)
            print(os.path.join(subdir, f))
        else:
            os.remove(os.path.join(subdir, f))
