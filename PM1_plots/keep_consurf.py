import os

rootdir = "/home/fryerv/sc_git_clone/snippets/PM1_plots/consurf_scores"

for subdir, dir, files in os.walk(rootdir):
    for f in files:
        if f == "consurf.grades" or f == "keep_consurf.py":
            print(f)
            print(os.path.join(subdir, f))
        else:
            os.remove(os.path.join(subdir, f))
