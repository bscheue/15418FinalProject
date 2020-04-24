import math
import random


def neighbors((x, y)):
  return {(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)}


def neighborhood(cluster):
  res = set()
  for point in cluster:
    res.update(neighbors(point))
  return res


def prefixsum(xs, start):
  accx = start[0]
  accy = start[1]
  res = [None for i in range(len(xs))]
  for i in range(len(xs)):
    accx += xs[i][0]
    accy += xs[i][1]
    res[i] = (accx, accy)
  return res


# from page 14
def chooseK(rc):
  factor = 3
  return factor * rc**2


# from page 14
def chooseW(M):
  factor = 3
  epsilon = .01
  return factor * M**(1 + epsilon)


def step1(rc, M):
  rb = rc + 2
  k = int(chooseK(rc))
  w = int(chooseW(M))
  return (rb, k, w)


def createStep(n):
  if n == 1:
    return (1, 0)
  elif n == 2:
    return (0, 1)
  elif n == 3:
    return (-1, 0)
  elif n == 4:
    return (0, -1)


def createWalk(rb, k, start):
  dirs = [random.randint(1, 4) for i in range(k)]
  steps = [createStep(n) for n in dirs]
  locations = prefixsum(steps, start)
  return locations


def createStart(rb):
  deg = random.randint(1, 360)
  return (rb * int(math.sin(deg)), rb * int(math.cos(deg)))


def step2(rb, k, w):
  # return walks
  return [createWalk(rb, k, createStart(rb)) for i in range(w)]


def step3(walks, cluster):
  res = [None for i in range(len(walks))]
  nbhd = neighborhood(cluster)
  for i in range(len(walks)):
    for j in range(len(walks[i])):
      if walks[i][j] in nbhd:
        res[i] = j
        break
  return res


def step4(res, walks):
  for i in range(len(walks)):
    for j in range(i):
      if res[j] is None:
        break
      i_walk = walks[i]
      for k in range(len(i_walk)):
        if walks[j][res[j]] == i_walk[k]:
          return i
  return None


def step5(cluster, res, walks, k):
  k = k if k is not None else len(res)
  for i in range(k):
    if res[i] is not None:
      cluster.add(walks[i][res[i]])


def doBatch(cluster):
  rc = 1
  M = 1
  (rb, k, w) = step1(rc, M)
  walks = step2(rb, k, w)
  res = step3(walks, cluster)
  k = step4(res, walks)
  step5(cluster, res, walks, k)
  print(cluster)


def doSimulation():
  cluster = set()
  cluster.add((0, 0))
  doBatch(cluster)


if __name__ == "__main__":
  doSimulation()
