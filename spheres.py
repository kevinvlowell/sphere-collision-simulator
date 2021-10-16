# Copyright 2021 Samir Farhat Dominguez safarhat@bu.edu
# Copyright 2021 Fahad Farid fahd@bu.edu
# Kevin, please add your info!
import sys

sphere_list = []
events = []

# Our main
def run_sim():
    for sphere in sphere_list:
        print(sphere)



def report_event(event):
    print("event report")



initial_conditions = sys.argv
if(len(initial_conditions) != 3):
    raise Exception("Please input initial conditions properly: radius duration")

universe_radius = float(initial_conditions[1])
sim_duration = float(initial_conditions[2])

print("Please enter the mass, radius, x/y/z position, x/y/z velocity and name of each sphere")
print("When complete, use EOF / Ctrl-D to stop entering")

# Take input
for line in sys.stdin:
    sphere_list.append(line)


print("Here are the initial conditions.")
print("universe radius " + str(universe_radius))
print("end simulation " + str(sim_duration))

run_sim()

print("energy: ")
print("momentum: ")
print("Here are the events")

for event in events:
    report_event(event)