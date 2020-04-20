import math
import random


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
  factor = 1
  return factor * rc**2


# from page 14
def chooseW(M):
  factor = 1
  epsilon = .01
  return factor * M**(1 + epsilon)


def step1(rc):
  rb = rc + 2
  k = chooseK(rc)
  w = chooseW(M)
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
  return [createWalk(rb, k, createStart(rb)) for i in range(w)]
