direc="/home/mitsuki/altorf/denovo/data"

def get_strain_lst(target):
    fp = "{}/{}/strain.lst".format(direc, target)
    return [s.strip() for s in open(fp, 'r').readlines()]

if __name__=="__main__":
    print(get_strain_lst("synechococcus"))
