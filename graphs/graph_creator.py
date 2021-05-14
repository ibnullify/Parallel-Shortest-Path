import sys
import random 


NUM_VERTICES = int(sys.argv[1])

OUT_FILE = "graph{}".format(NUM_VERTICES)


def make_graph():
    with open(OUT_FILE, 'w') as f:
        f.write(str(NUM_VERTICES) + "\n")
        for i in range(NUM_VERTICES-1):
            for j in range(i+1):
                cost = random.randint(0,200)
                f.write(str(cost) + ' ')
            if i != NUM_VERTICES-2: f.write("\n")
    f.close()



make_graph()