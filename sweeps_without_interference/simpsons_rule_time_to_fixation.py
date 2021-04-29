#This is to write up Simpson'r rule of numerical integration for any function f(x) continuous in interval [a,b]:
#Choose an even integer n, which will be the number of intervals we divide b to a in, given by x0, x1, x2 -- xn
#d = (b-a)/n
#Then, integral (I) of f(x) from x=a to x=b will be:
#I = d/3*[f(x0) + 4f(x1) + 2f(x2) + .. + 4f(xn-1) + f(xn)]
#The pattern is [1, 4, 2, 4, .., 2,4,...2, 4, 1]
#Example of how to run:
#python simpsons_rule.py -numIntervals 500 -Ne 10000 -Gamma 10

import sys
import math
import argparse

#parsing user given constants
parser = argparse.ArgumentParser(description='Population parameters etc')
parser.add_argument('-numIntervals', dest = 'numIntervals', action='store', nargs = 1, type = int, default=500, help = 'number of intervals, should be even. Larger this is, more accurate the answer.')
parser.add_argument('-Ne', dest = 'Ne', action='store', nargs = 1, type = int, default=10000, help = 'Wright-Fisher effective population size')
parser.add_argument('-Gamma', dest = 'Gamma', action='store', nargs = 1, type = float, default=10, help = 'gamma = 2Nes')
args = parser.parse_args()
n = int(args.numIntervals[0])
N = int(args.Ne[0])
gamma = int(args.Gamma[0])
#n = int(sys.argv[1])#number of intervals, should be even

if n > 0 and n%2==0:
    print ("Integration over " + str(n) + " intervals")
else:
    print ("n needs to be even and greater than zero")

#define constants
#gamma = 9.5 #50.0
S = gamma/2.0
#N = 10000 #1000000.0
s = float(S)/float(N)

def J1(x):#needs to be integrated from p to 1
    N = 2.0*(math.exp(2.0*S*x)-1.0)*(math.exp(-2.0*S*x)-math.exp(-2.0*S))
    D = s*(1.0-math.exp(-2.0*S))*x*(1.0-x)
    y = float(N)/float(D)
    return(y)

def J2(x):
    N = 2.0*(math.exp(2.0*S*x)-1.0)*(1.0-math.exp(-2.0*S*x))
    D = s*(1.0-math.exp(-2.0*S))*x*(1.0-x)
    y = float(N)/float(D)
    return(y)

def U(x):
    N = 1 - math.exp(-2*S*x)
    D = 1 - math.exp(-2*S)
    y = float(N)/float(D)
    return(y)

def J2_limit(x):
    N = 2.0**2.0*S*(1-math.exp(-2.0*S*x))
    D = s*(1.0-math.exp(-2.0*S))*(1.0-x)
    y = N/float(D)
    return(y)

def time_to_fixation(x):
    u_x = U(x)
    t = J1(x) + (((1-u_x)/u_x)*J2(x))
    return(t)

def get_sequence(num_intervals):
    l_seq = [1,4]
    while len(l_seq) <= (num_intervals+1)-3:
        l_seq.append(2)
        l_seq.append(4)
    l_seq.append(1)
    return(l_seq)

#Integrate J1 from p1 to p2:
p1 = 1.0/(2*N)
p2 = float((2*N)-1)/float(2*N)
d = (p2 - p1)/float(n)
l_x = [p1]
i = 1
while i <= n:
    l_x.append(p1 + (float(i)*d))
    i = i+1
print(l_x)
l_coeff = get_sequence(n)
print (l_coeff)

#perform a check here:
if len(l_x) != len(l_coeff):
    print ("something is wrong.")
else:
    print ("length of vectors is: " + str(len(l_x)))
#Get J1:
I1 = 0.0
j = 0
while j < n+1:
    I1 = I1 + ((d/3.0)*float(l_coeff[j])*J1(l_x[j]))
    j = j + 1
print (I1)

#Integrate J2 from 0 to p1, in the limit of p1 tending to zero:
d = (p1 - 0.0)/float(n)
l_x = [0.0]
i = 1
while i <= n:
    l_x.append(0.0 + (float(i)*d))
    i = i + 1
#Get J2:
I2 = 0.0
j = 0
while j < n+1:
    I2 = I2 + ((d/3.0)*float(l_coeff[j])*J2_limit(l_x[j]))
    j = j + 1
print (I2)
#Get time to fixation:
u_p = U(p1)
t = I1 + (((1-u_p)/u_p)*I2)
print("number of generations: " + str(t))
print("number of 2N generations: " + str(t/(2.0*N)))
print ("done")

