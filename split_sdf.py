from sys import argv

def split_SDF(filename):
    with open(filename,'r') as f:
        tmp = []
        index = 0
        for line in f:
            if line != "$$$$\n":
                tmp.append(line)
            else:
                tmp.append(line)
                with open('{0}_{1}'.format(index,filename),'w+') as wt:
                    for ds in tmp:
                        wt.write(ds)
                index += 1
                tmp.clear()

if __name__ == "__main__":
    split_SDF(argv[1])