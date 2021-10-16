# Copyright 2021 Samir Farhat Dominguez safarhat@bu.edu
# Copyright 2021 Fahad Farid fahd@bu.edu
# Copyright 2021 Kevin Vogt-Lowell kjv@bu.edu

import sys
import math

class Sphere:

    def __init__(self, attributes_string):

        attributes_list = attributes_string.split()

        if len(attributes_list) == 9:
            self.mass = float(attributes_list[0])
            self.radius = float(attributes_list[1])

            self.pos_x = float(attributes_list[2])
            self.pos_y = float(attributes_list[3])
            self.pos_z = float(attributes_list[4])

            self.vel_x = float(attributes_list[5])
            self.vel_y = float(attributes_list[6])
            self.vel_z = float(attributes_list[7])

            self.name = attributes_list[8]

        else:
            raise Exception("Attributes not properly specified for sphere creation. Please check input and try again.")



sphere_strings = []
sphere_list = []
events = []

# Our main
def run_sim(radius, duration):
    for string in sphere_strings:
        sphere_list.append(Sphere(string))

    energy = 0
    momentum = [0, 0, 0] #list entries correspond to x, y, z directions, respectively

    #print out all initial conditions
    print("Here are the initial conditions.")
    print("universe radius " + str(radius))
    print("end simulation " + str(duration))
    for sphere in sphere_list:
        print(sphere.name, "m="+str(sphere.mass), "R="+str(sphere.radius),\
        "p=(" + str(sphere.pos_x) + "," + str(sphere.pos_y) + "," + str(sphere.pos_z) + ")",\
        "v=(" + str(sphere.vel_x) + "," + str(sphere.vel_y) + "," + str(sphere.vel_z) + ")")

        #calculate initial energy (0.5mv^2, summed over spheres) and momentum (mv vector)
        velocity_magnitude = (sphere.vel_x**2 + sphere.vel_y**2 + sphere.vel_z**2)
        energy += 0.5 * sphere.mass * velocity_magnitude**2

        momentum[0] += sphere.mass * sphere.vel_x
        momentum[1] += sphere.mass * sphere.vel_y
        momentum[2] += sphere.mass * sphere.vel_z

    print("energy:", str(energy))
    print("momentum: (" + str(momentum[0]) + "," + str(momentum[1]) + "," + str(momentum[2]) + ")")
    print("\nHere are the events:")

    for event in events:
        report_event(event)


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

run_sim(universe_radius, sim_duration)
