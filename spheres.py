# Copyright 2021 Samir Farhat Dominguez safarhat@bu.edu
# Copyright 2021 Fahad Farid fahd@bu.edu
# Copyright 2021 Kevin Vogt-Lowell kjv@bu.edu

import sys

class Sphere:

    def __init__(self, attributes_string):

        attributes_list = attributes_string.split()

        if len(attributes_list) == 9:
            self.mass = attributes_list[0]
            self.radius = attributes_list[1]

            self.pos_x = attributes_list[2]
            self.pos_y = attributes_list[3]
            self.pos_z = attributes_list[4]

            self.vel_x = attributes_list[5]
            self.vel_y = attributes_list[6]
            self.vel_z = attributes_list[7]

            self.name = attributes_list[8]

        else:
            raise Exception("Attributes not properly specified for sphere creation. Please check input and try again.")



sphere_strings = []
sphere_list = []
events = []

# Our main
def run_sim():
    for string in sphere_strings:
        sphere_list.append(Sphere(string))

    print("Sphere list: ", sphere_list)


def report_event(event):
    print("event report")



initial_conditions = sys.argv
if(len(initial_conditions) != 3):
    raise Exception("Please input initial conditions properly: radius, duration")

universe_radius = float(initial_conditions[1])
sim_duration = float(initial_conditions[2])

print("Please enter the mass, radius, x/y/z position, x/y/z velocity and name of each sphere")
print("When complete, use EOF / Ctrl-D to stop entering")

# Take input
for line in sys.stdin:
    sphere_strings.append(line)


print("Here are the initial conditions:")
print("Universe radius: " + str(universe_radius))
print("Simulation duration: " + str(sim_duration))

run_sim()

print("energy: ")
print("momentum: ")
print("Here are the events:")

for event in events:
    report_event(event)
