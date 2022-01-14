from argparse import ArgumentParser
import time
import numpy
import os
import psutil

DELTA = 30

ALPHA = {
    'A': {
        'A': 0,
        'C': 110,
        'G': 48,
        'T': 94,
    },
    'C': {
        'A': 110,
        'C': 0,
        'G': 118,
        'T': 48,
    },
    'G': {
        'A': 48,
        'C': 118,
        'G': 0,
        'T': 110,
    },
    'T': {
        'A': 94,
        'C': 48,
        'G': 110,
        'T': 0,
    },
}


def sequence_alignment(s1, s2):
    DP = numpy.zeros((len(s1) + 1, len(s2) + 1))

    for i in range(1, len(s1) + 1):
        DP[i][0] = DELTA * i

    for j in range(1, len(s2) + 1):
        DP[0][j] = DELTA * j

    s1 = 'x' + s1
    s2 = 'x' + s2

    for i in range(1, len(s1)):
        for j in range(1, len(s2)):
            DP[i][j] = min(
                DP[i][j - 1] + DELTA,
                DP[i - 1][j] + DELTA,
                DP[i - 1][j - 1] + ALPHA[s1[i]][s2[j]],
            )

    return DP


def backtrack(s1, s2, DP):
    s1 = 'x' + s1
    s2 = 'x' + s2

    i = len(s1) - 1
    j = len(s2) - 1

    t1 = []
    t2 = []

    while i > 0 and j > 0:
        x = DP[i][j - 1] + DELTA
        y = DP[i - 1][j] + DELTA
        z = DP[i - 1][j - 1] + ALPHA[s1[i]][s2[j]]

        m = min(x, y, z)

        if z == m:
            t1.append(s1[i])
            t2.append(s2[j])
            i -= 1
            j -= 1
        elif y == m:
            t1.append(s1[i])
            t2.append('_')
            i -= 1
        else:
            t1.append('_')
            t2.append(s2[j])
            j -= 1

    while i > 0:
        t1.append(s1[i])
        t2.append('_')
        i -= 1

    while j > 0:
        t1.append('_')
        t2.append(s2[j])
        j -= 1

    t1 = ''.join(reversed(t1))
    t2 = ''.join(reversed(t2))

    return t1, t2


def read_input(filename):
    with open(filename, 'r') as f:
        data = f.readlines()

    s1 = data[0].strip()
    js = []

    for i, dat in enumerate(data[1:]):
        try:
            js.append(int(dat.strip()))
        except ValueError:
            break

    s2 = data[i + 1].strip()
    ks = []

    for dat in data[i + 2:]:
        ks.append(int(dat.strip()))

    for j in js:
        s1 = s1[:j + 1] + s1 + s1[j + 1:]

    for k in ks:
        s2 = s2[:k + 1] + s2 + s2[k + 1:]

    return s1, s2


def write_output(t, m, t1, t2, cost):
    with open('output.txt', 'w') as f:
        f.write('{} {}\n'.format(t1[:50], t1[-50:]))
        f.write('{} {}\n'.format(t2[:50], t2[-50:]))
        f.write('{:.3f}\n'.format(t))
        f.write('{}'.format(m))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_file')

    args = parser.parse_args()

    s1, s2 = read_input(args.input_file)


    start_time = time.time()

    DP = sequence_alignment(s1, s2)
    t1, t2 = backtrack(s1, s2, DP)

    total_time = time.time() - start_time
    total_memory = psutil.Process(os.getpid()).memory_info().rss / 1024


    write_output(total_time, total_memory, t1, t2, DP[-1][-1])
