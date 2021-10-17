# Copyright 2021 Samir Farhat Dominguez safarhat@bu.edu
# Copyright 2021 Fahad Farid fahd@bu.edu
# Copyright 2021 Kevin Vogt-Lowell kjv@bu.edu

import sys, math

######################################################################################################
######################################################################################################
######################################################################################################


# Sphere Class
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


# Event Class
class Event:

    #assumes we would pass the list of spheres after adjusting the velocities in the simulation code according to the event
    #would also calculate new energy and momentum for arguments immediately prior to creating event
    def __init__(self, event_time, event_type, sphere_a, sphere_b, all_spheres, energy, momentum):

        self.sphere_names_involved = [sphere_a.name]
        self.new_sphere_statuses = all_spheres #can (hopefully) index this list to get saved attributes of all spheres at the time of the event

        self.timestamp = float(event_time)
        self.type = event_type

        if self.type == "colliding":
            self.sphere_names_involved.append(sphere_b.name)

        self.energy = energy
        self.momentum = momentum


######################################################################################################
######################################################################################################
######################################################################################################


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


def compute_collision_time(s1, s2):
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

		return minimum_pos_quadratic_soln(a,b,c)
	return -2


def compute_reflection_time(sphere,radius):
	#check when each sphere would hit the wall if moving and not blocked
	##checking to see if sphere is moving
	if not(sphere.vel_x == 0 and sphere.vel_y == 0 and sphere.vel_z == 0):
		a = sphere.vel_x**2 + sphere.vel_y**2 + sphere.vel_z**2
		b = (2*sphere.pos_x)*sphere.vel_x + (2*sphere.pos_y)*sphere.vel_y + (2*sphere.pos_z)*sphere.vel_z
		c = sphere.pos_x**2 + sphere.pos_y**2 + sphere.pos_z**2 - (radius - sphere.radius)**2
		return minimum_pos_quadratic_soln(a,b,c)
	return -2


# get collsion results
def compute_collision(sphere_a, sphere_b):
	v1_x = sphere_a.vel_x
	v1_y = sphere_a.vel_y
	v1_z = sphere_a.vel_z
	v1_x = sphere_b.vel_x
	v1_y = sphere_b.vel_y
	v1_z = sphere_b.vel_z

	v1_x_prime = 0

	pass


def compute_reflection(sphere):
	pass


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

    # begin simulation
    ## determine the time of the first event by evaluating all possible collision times from starting
    # positions, and choosing event that happens soonest
    nearest_event_time = -1
    next_colliding_pair_indices = []

    # check sphere-to-sphere collision times
    for i in range(len(sphere_list)-1):
        s1 = sphere_list[i]
        j = i+1

        while j < len(sphere_list):
            s2 = sphere_list[j]
            current_event_time = compute_collision_time(s1,s2)
            # compare current event time to value of nearest_event_time to see if this event would happen sooner, and update nearest_event_time if so
            if (current_event_time < nearest_event_time or nearest_event_time == -1) and current_event_time != -2:
                nearest_event_time = current_event_time
                next_event_type = "colliding"
                next_colliding_pair_indices = [i, j]
            # increment j to look at next possible sphere pair
            j += 1

    # check sphere-to-wall collision times - WILL NEED TESTING ONCE VELOCITIES ARE BEING ADJUSTED;
    # Notes (3.4)on assignment says objects are never initially in a overlapping or colliding state
    for sphere in sphere_list:
            # check when each sphere would hit the wall if moving and not blocked
            ## checking to see if sphere is moving
            if not(sphere.vel_x == 0 and sphere.vel_y == 0 and sphere.vel_z == 0):
                current_event_time = compute_reflection_time(sphere,radius)
                # compare current event time to value of nearest_event_time to see if this event would happen sooner, and update nearest_event_time if so
                if (current_event_time < nearest_event_time or nearest_event_time == -1) and current_event_time != -2:
                    nearest_event_time = current_event_time
                    next_event_type = "reflecting"
                    next_reflecting_sphere = sphere


    #make adjustments to sphere velocities based on next occurring event (use nearest_event_time and next_event_type), and add event to event list using event class
    print("\nnearest_event_time:", nearest_event_time)
    # end simulation once nearest event is past the existence of the universe
    if(nearest_event_time > sim_duration):
         sys.exit(0)


    if next_event_type == "reflecting":
        # Compute new velocities, system energy, and system momentum


        #create an event and append to event list
        current_event = Event(nearest_event_time, next_event_type, next_reflecting_sphere, None, sphere_list, energy, momentum)
        event_list.append(current_event)

    elif next_event_type == "colliding":
        # Compute new velocities, system energy, and system momentum
        colliding_sphere1 = sphere_list[next_colliding_pair_indices[0]]
        colliding_sphere2 = sphere_list[next_colliding_pair_indices[1]]



        # create an event and append to event list
        current_event = Event(nearest_event_time, next_event_type, colliding_sphere1, colliding_sphere2, sphere_list, energy, momentum)
        event_list.append(current_event)

    else:
        pass


    print("\nHere are the events:")

    for event in event_list:
        report_event(event)


def report_event(event):
    print("event report")


######################################################################################################
######################################################################################################
######################################################################################################


# Script starts here
sphere_strings = []
sphere_list = []
event_list = []


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
