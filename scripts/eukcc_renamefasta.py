# quick script to make sure names in fasta files are not stupid for gmes

import sys

inf = sys.argv[1]
outf = sys.argv[2]

nms = []
with open(outf, "w") as o:
    with open(inf) as f:
        for line in f:
            if line.startswith(">"):
                l = line.split()
                N = l[0].strip()
                n = N
                i = 0
                while n in nms:
                    n = "{}_{}".format(N ,i)
                    i += 1
                o.write("{}\n".format(n))
                nms.append(n)
            else:
                o.write(line)