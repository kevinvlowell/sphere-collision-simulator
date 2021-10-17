"""This program generates a collection
of results (saved in text files) for the spheres problem by
repeatedly running spheres.py with the 
parameters specified in the list "spheres_tests"
defined below.

The file name is partly based on "myid"

My results were created with "myid" set to jbc

You should change this variable so you can
compare your results with mine.

"""

myid = "jbc"

from subprocess import run, PIPE
import os

try:
  os.mkdir('examples')
except:
  pass
finally:
  print('results will be in examples/')

spheres_tests=[
('tennisball',100,100),
('tennisball',10.41,200),
('zmover',100,300),
('xyzmovers',200,200),
('three',120,200), 
('threeonxaxis',50,400),
('diagonal',100,200),
('randompool',200,500),
('lineup',2000,200),
('smashup',2000,500),
('slow',2000,1e84),
('glancing',40,200)]


for basename,radius,N in spheres_tests:
    with open(f'examples/{basename}.txt') as g:
        text = g.read()

    T = run(['python','spheres.py',str(radius),str(N)],input=text,stdout=PIPE,
          universal_newlines=True)
    
    with open(f"examples/{basename}_{radius}_{N}_{myid}.txt",'w') as out:
        out.write(T.stdout)
    
