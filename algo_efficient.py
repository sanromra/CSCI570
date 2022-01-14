import time
from argparse import ArgumentParser
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


def sequence_alignment_efficient(s1, s2):

    s1 = 'x' + s1
    s2 = 'x' + s2

    cols = numpy.zeros((len(s1), 2))

    for i in range(len(s1)):
        cols[i][0] = i * DELTA

    cols[0][1] = DELTA

    j = 1

    while j < len(s2):
        for i in range(1, len(s1)):
            cols[i][1] = min(
                cols[i-1][1] + DELTA,
                cols[i][0] + DELTA,
                cols[i-1][0] + ALPHA[s1[i]][s2[j]]
            )

        j += 1

        cols[:, 0] = cols[:, 1]
        cols[:, 1] = 0
        cols[0][1] = j*DELTA

    return cols


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


def divide_and_conquer(s1, s2):
    m = len(s1)
    n = len(s2)

    if m <= 2 or n <= 2:
        DP = sequence_alignment(s1, s2)
        return backtrack(s1, s2, DP)

    min_idx = numpy.inf
    min_sum = numpy.inf

    c1 = sequence_alignment_efficient(s1, s2[:n // 2])
    c2 = sequence_alignment_efficient(s1[::-1], s2[n // 2:][::-1])

    c2 = c2[::-1]

    for idx in range(c2.shape[0]):
        s = c1[idx][0] + c2[idx][0]
        if s < min_sum:
            min_sum = s
            min_idx = idx

    l1, l2 = divide_and_conquer(s1[:min_idx], s2[:n // 2])
    r1, r2 = divide_and_conquer(s1[min_idx:], s2[n // 2:])

    return l1 + r1, l2 + r2


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


def write_output(t, m, t1, t2):
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

    e1, e2 = divide_and_conquer(s1, s2)

    total_time = time.time() - start_time
    total_memory = psutil.Process(os.getpid()).memory_info().rss / 1024

    write_output(total_time, total_memory, e1, e2)
