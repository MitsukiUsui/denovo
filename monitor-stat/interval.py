#!/usr/bin/env python3

class Interval:
    def __init__(self, start, end):
        assert start <= end
        self.start = start
        self.end = end
        
    def __str__(self):
        return "[{},{})".format(self.start, self.end)
    
def interval_justsum(interval_lst):
    ret = 0
    for interval in interval_lst:
        ret += interval.end - interval.start
    return ret

def interval_sum(interval_lst):
    pos_lst = []
    for interval in interval_lst:
        pos_lst.append((interval.start, 0))
        pos_lst.append((interval.end, 1))
    pos_lst=sorted(pos_lst, key=lambda x: (x[0], x[1]))
    
    prv, cover, count = 0, 0, 0
    for pos in pos_lst:
        move = pos[0] - prv
        if count > 0:
            cover += move

        if pos[1] == 0:
            count += 1
        elif pos[1] == 1:
            count -= 1
        else:
            assert False
        
        assert count >= 0
        prv = pos[0]
    
    return cover

def interval_not(interval_lst, length):
    pos_lst = []
    for interval in interval_lst:
        pos_lst.append((interval.start, 0))
        pos_lst.append((interval.end, 1))
    pos_lst=sorted(pos_lst, key=lambda x: (x[0], x[1]))
    
    ret_lst = []
    prv, count = 0, 0
    for pos in pos_lst:
        if count == 0:
            assert pos[1] == 0
            ret_lst.append(Interval(prv, pos[0]))

        if pos[1] == 0:
            count += 1
        elif pos[1] == 1:
            count -= 1
        else:
            assert False
        assert count >= 0
        prv = pos[0]
    
    if prv < length: #when CDS overlap breaks in genome, gff.end is larger then its genome length
        ret_lst.append(Interval(prv, length))
    return ret_lst

def interval_and(intervals1, intervals2):
    pos_lst = []
    for interval in intervals1:
        pos_lst.append((interval.start, 0, 0))
        pos_lst.append((interval.end, 1, 0))
    for interval in intervals2:
        pos_lst.append((interval.start, 0, 1))
        pos_lst.append((interval.end, 1, 1))
    pos_lst=sorted(pos_lst, key=lambda x: (x[0], x[1]))
    
    prv, cover =  0, 0
    count = [0, 0]
    
    for pos in pos_lst:
        move = pos[0] - prv
        if min(count) > 0:
            cover += move

        if pos[1] == 0:
            count[pos[2]] += 1
        elif pos[1] == 1:
            count[pos[2]] -= 1        
        assert min(count) >= 0
        prv = pos[0]
    
    return cover
