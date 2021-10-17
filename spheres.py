# Copyright 2021 Samir Farhat Dominguez safarhat@bu.edu
# Copyright 2021 Fahad Farid fahd@bu.edu
# Copyright 2021 Kevin Vogt-Lowell kjv@bu.edu

import sys, math

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



# EVENT CLASS
class Event:

    def __init__(self,log,sphere_a,sphere_b,spheres):
        taken_spheres = {sphere_a.name}
        log_other_spheres = []

        self.timestamp = float(log[0])
        self.type = log[1]

        if(self.type == "colliding"):
            self.colliding = {log[2], log[3]}
            self.colliding_1 = sphere_a
            self.colliding_2 = sphere_b
            taken_spheres.append(sphere_b.name)
        elif(self.type == "reflecting"):
            self.colliding_1 = sphere_a

        for sphere in spheres:
            if sphere.name not in taken_spheres:
                log_other_spheres.append(sphere)

        self.energy = log[4]
        self.momentum  = log[5]




def minimum_pos_quadratic_soln(a, b, c):
    #returns the smallest positive value that solves the quadratic equation

    t1 = (-b + math.sqrt(b**2 - 4*a*c))/(2*a)
    t2 = (-b - math.sqrt(b**2 - 4*a*c))/(2*a)

    #print("T1 =",t1)
    #print("T2 =",t2)

    if t1 >= 0 and t2 >= 0:
        return min(t1,t2)
    elif t1 >= 0 and t2 < 0:
        return t1
    elif t1 < 0 and t2 >= 0:
        return t2
    else:
        raise Exception("Error: no positive collision times.")




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
    print("\nHere are the initial conditions.")
    print("universe radius " + str(radius))
    print("end simulation " + str(int(duration)))
    
    for sphere in sphere_list:
        print(sphere.name, "m="+str(sphere.mass), "R="+str(sphere.radius),\
        "p=(" + str(int(sphere.pos_x)) + "," + str(int(sphere.pos_y)) + "," + str(int(sphere.pos_z)) + ")",\
        "v=(" + str(int(sphere.vel_x)) + "," + str(int(sphere.vel_y)) + "," + str(int(sphere.vel_z)) + ")")

        #calculate initial energy (0.5mv^2, summed over spheres) and momentum (mv vector)
        velocity_magnitude = (sphere.vel_x**2 + sphere.vel_y**2 + sphere.vel_z**2)
        energy += 0.5 * sphere.mass * velocity_magnitude**2

        momentum[0] += sphere.mass * sphere.vel_x
        momentum[1] += sphere.mass * sphere.vel_y
        momentum[2] += sphere.mass * sphere.vel_z

    print("energy:", str(int(energy)))
    print("momentum: (" + str(int(momentum[0])) + "," + str(int(momentum[1])) + "," + str(int(momentum[2])) + ")")

    #begin simulation
    ##determine the time of the first event by evaluating all possible collision times from starting 
    # positions, and choosing event that happens soonest
    nearest_event_time = -1
    next_colliding_pair_indices = []

    #check sphere-to-sphere collision times
    for i in range(len(sphere_list)-1):

        s1 = sphere_list[i]
        j = i+1

        while j < len(sphere_list):

            s2 = sphere_list[j]

            #only consider sphere pairs that are approaching each other, aka if <(r1-r2),(v1-v2)> less than 0
            pos_diff_x = s1.pos_x-s2.pos_x
            pos_diff_y = s1.pos_y-s2.pos_y
            pos_diff_z = s1.pos_z-s2.pos_z

            vel_diff_x = s1.vel_x-s2.vel_x
            vel_diff_y = s1.vel_y-s2.vel_y
            vel_diff_z = s1.vel_z-s2.vel_z

            if (pos_diff_x*vel_diff_x + pos_diff_y*vel_diff_y + pos_diff_z*vel_diff_z) < 0:

                #calculate quadratic equation coefficients based on current sphere pair

                ## a = <v1,v1> - 2<v1,v2> + <v2,v2>
                a = (s1.vel_x**2 + s1.vel_y**2 + s1.vel_z**2) - 2*(s1.vel_x*s2.vel_x + s1.vel_y*s2.vel_y + s1.vel_z*s2.vel_z) + (s2.vel_x**2 + s2.vel_y**2 + s2.vel_z**2)

                ## b = <p1,v1> - <p1,v2> - <p2,v1> + <p2,v2>
                b = (s1.pos_x*s1.vel_x + s1.pos_y*s1.vel_y + s1.pos_z*s1.vel_z) - (s1.pos_x*s2.vel_x + s1.pos_y*s2.vel_y + s2.pos_z*s1.vel_z) - (s2.pos_x*s1.vel_x + s2.pos_y*s1.vel_y + s2.pos_z*s1.vel_z) + (s2.pos_x*s2.vel_x + s2.pos_y*s2.vel_y + s2.pos_z*s2.vel_z)

                ## c = <p1,p1> - 2<p1,p2> + <p2,p2> - (R1 + R2)**2
                c = (s1.pos_x**2 + s1.pos_y**2 + s1.pos_z**2) - 2*(s1.pos_x*s2.pos_x + s1.pos_y*s2.pos_y + s1.pos_z*s2.pos_z) + (s2.pos_x**2 + s2.pos_y**2 + s2.pos_z**2) - (s1.radius + s2.radius)**2

                # solve quadratic equation to find collision times for current sphere pair
                current_event_time = minimum_pos_quadratic_soln(a,b,c)

                #compare current event time to value of nearest_event_time to see if this event would happen sooner, and update nearest_event_time if so
                if current_event_time < nearest_event_time or nearest_event_time == -1:
                    nearest_event_time = current_event_time
                    next_event_type = "colliding"
                    next_colliding_pair_indices = [i, j]

            #increment j to look at next possible sphere pair
            j += 1

    #check sphere-to-wall collision times - WILL NEED TESTING ONCE VELOCITIES ARE BEING ADJUSTED; 
    # Notes (3.4)on assignment says objects are never initially in a overlapping or colliding state
    for sphere in sphere_list:

            #check when each sphere would hit the wall if moving and not blocked

            ##checking to see if sphere is moving
            if not(sphere.vel_x == 0 and sphere.vel_y == 0 and sphere.vel_z == 0):

                a = sphere.vel_x**2 + sphere.vel_y**2 + sphere.vel_z**2

                b = (2*sphere.pos_x)*sphere.vel_x + (2*sphere.pos_y)*sphere.vel_y + (2*sphere.pos_z)*sphere.vel_z

                c = sphere.pos_x**2 + sphere.pos_y**2 + sphere.pos_z**2 - (radius - sphere.radius)**2

                current_event_time = minimum_pos_quadratic_soln(a,b,c)

                #compare current event time to value of nearest_event_time to see if this event would happen sooner, and update nearest_event_time if so
                if current_event_time < nearest_event_time or nearest_event_time == -1:
                    nearest_event_time = current_event_time
                    next_event_type = "reflecting"
                    next_reflecting_sphere = sphere


    #make adjustments to sphere velocities based on next occurring event (use nearest_event_time and next_event_type), and add event to event list using event class
    print("\nnearest_event_time:", nearest_event_time)

    if next_event_type == "reflecting":
        pass
    elif next_event_type == "colliding":
        pass
    else:
        pass


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
