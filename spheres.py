# Copyright 2021 Samir Farhat Dominguez safarhat@bu.edu
# Copyright 2021 Fahad Farid fahd@bu.edu
# Kevin, please add your info!

import sys

# Our main
def run_sim(params):
    print(params)


initial_conditions = sys.argv
if(len(initial_conditions) != 2):
    print("Please input initial conditions properly: radius duration")

try:
    universe_radius = int(initial_conditions[0])
except:
    print("Input universe radius correctly please")

try:
    sim_duration = int(initial_conditions[1])
except:
    print("Input sim duration correctly please")


# Take input
for line in sys.stdin:
    run_sim(line)
