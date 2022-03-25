import numpy as np
import argparse
import sys

def readFASTA(filename):
    fullString = ''
    with open(filename) as file:
        lines = file.readlines()
        for i in range(1, len(lines)):
            tempLine = ''.join(lines[i].split())
            fullString += tempLine
    return fullString

def readBLOSUM(filename):
    outerDict = {}
    with open(filename) as f:
        lines = f.readlines()
        acids = lines[0].split()
        for i in range(1, len(lines)):
            splitLine = lines[i].split()
            dict = {}
            for j in range(1, len(splitLine)):
                dict[acids[j-1]] = int(splitLine[j])
            outerDict[splitLine[0]] = dict
    return outerDict

def makeMatrix(blosum, v, w, penalty, penalty_ext):
    qryMatrix = np.zeros(shape=(len(v), len(w)))
    sbjMatrix = np.zeros(shape=(len(v), len(w)))
    scoreMatrix = np.zeros(shape=(len(v)+1, len(w)+1))
    backtrack = np.zeros(shape=(len(v), len(w)))

    for i in range(1, len(v)):
        for j in range(1, len(w)):

            qryMatrix[i][j] = max(0, scoreMatrix[i-1][j] - penalty, qryMatrix[i-1][j] - penalty_ext)

            sbjMatrix[i][j] = max(0, scoreMatrix[i][j-1] - penalty, sbjMatrix[i][j-1] - penalty_ext)

            diagonal = scoreMatrix[i - 1][j - 1] + blosum[v[i]][w[j]]
            delete = qryMatrix[i - 1][j - 1] + blosum[v[i]][w[j]]
            insert = sbjMatrix[i - 1][j - 1] + blosum[v[i]][w[j]]

            scoreMatrix[i][j] = max(diagonal, 0, insert, delete)

            if scoreMatrix[i][j] == diagonal:
                backtrack[i][j] = 2
            elif scoreMatrix[i][j] == insert:
                backtrack[i][j] = 1
            elif scoreMatrix[i][j] == delete:
                backtrack[i][j] = 0
            elif scoreMatrix[i][j] == 0:
                backtrack[i][j] = 3

    return scoreMatrix, backtrack

def OutputLCS(backtrack, v, w, i, j, v_, w_):
    if i == 0 or j == 0:
        return v[i]+v_, w[j]+w_
    elif backtrack[i][j] == 3:
        return v_, w_
    elif backtrack[i][j] == 0:
        return OutputLCS(backtrack, v, w, i - 1, j, v[i]+v_, '-'+w_)
    elif backtrack[i][j] == 1:
        return OutputLCS(backtrack, v, w, i, j - 1, '-'+v_, w[j]+w_)
    else:
        return OutputLCS(backtrack, v, w, i - 1, j - 1, v[i]+v_, w[j]+w_)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare two protein sequences and find alignment')
    parser.add_argument('-p', metavar='--penalty', type=int, default=12, help='The initial penalty of gap opening')
    parser.add_argument('-pe', metavar='--penalty_ext', type=int, default=1, help='The penalty of gap extension')
    parser.add_argument('-f1', metavar='--file1', required=True, help='The name of the Query sequence file --REQUIRED')
    parser.add_argument('-f2', metavar='--file2', required=True, help='The name of the Subject sequence file --REQUIRED')

    p = parser.parse_args()
    penalty = p.p
    penalty_ext = p.pe
    file1 = p.f1
    file2 = p.f2
    v = readFASTA(file1)
    w = readFASTA(file2)

    blosum = readBLOSUM('BLOSUM62.txt')
    scoreMatrix, backtrack = makeMatrix(blosum, v, w, penalty, 1)
    score = np.amax(scoreMatrix)
    max_indexes = np.unravel_index(np.argmax(scoreMatrix), scoreMatrix.shape)
    v_, w_ = OutputLCS(backtrack, v, w, max_indexes[0], max_indexes[1], '', '')
    width = 60
    loops = int(len(v_)/width)
    print(score)
    for i in range(loops):
        q = 'Query: '
        s = 'Sbjct: '
        m = 'Match: '
        for j in range(width):
            q += v_[(i*width)+j]
        for j in range(width):
            if v_[(i*width)+j] == w_[(i*width)+j]:
                m += v_[(i*width)+j]
            elif v_[(i*width)+j] == '-' or w_[(i*width)+j] == '-':
                m += ' '
            elif blosum[v_[(i*width)+j]][w_[(i*width)+j]] > 0:
                m += '+'
            else:
                m += " "
        for j in range(width):
            s += w_[(i*width)+j]
        print(q)
        print(m)
        print(s)
        print()
    q = 'Query: '
    s = 'Sbjct: '
    m = 'Match: '

    for j in range((loops*width), len(v_)):
        q += v_[j]
    for j in range((loops*width), len(v_)):
        if v_[j] == w_[j]:
            m += v_[j]
        elif v_[j] == '-' or w_[j] == '-':
            m += ' '
        elif blosum[v_[j]][w_[j]] > 0:
            m += '+'
        else:
            m += " "
    for j in range((loops*width), len(v_)):
        s += w_[j]
    print(q)
    print(m)
    print(s)
    print()




