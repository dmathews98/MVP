

import numpy as np

def pbc(i):
    # Boolean so mod is only done on extraneous points needing image checked
    if (i > int(4-1)) or (i < 0):
        image_i = i%4

        return int(image_i)

    else:
        return i

def main():
    lattice = np.array([[1, 2, 3, 4],
                        [5, 6, 7, 8],
                        [9, 10, 11, 12],
                        [13, 14, 15, 16]])

    print(lattice)

    for i in range(4):
        for j in range(4):
            print(lattice[i][j])
            print('\n')
            for n in range(-1, 2):
                for m in range(-1, 2):
                    if n == 0 and m == 0:
                        continue
                    print(lattice[pbc(i+n)][pbc(j+m)])
            print('\n\n')

main()
