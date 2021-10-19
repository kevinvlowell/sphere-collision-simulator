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

    def __init__(self, event_time, event_type, sphere_a, sphere_b, all_spheres, energy, momentum):

        self.sphere_names_involved = [sphere_a.name]
        self.sphere_status_logs = []

        self.timestamp = float(event_time)
        self.type = event_type

        if self.type == "colliding":
            self.sphere_names_involved.append(sphere_b.name)

        for sphere in all_spheres:
            self.sphere_status_logs.append(f"{sphere.name} m={sphere.mass:g} R={sphere.radius:g} p=({sphere.pos_x:g},{sphere.pos_y:g},{sphere.pos_z:g}) v=({sphere.vel_x:g},{sphere.vel_y:g},{sphere.vel_z:g})")

        self.energy = energy
        self.momentum = momentum


######################################################################################################
######################################################################################################
######################################################################################################


def minimum_pos_quadratic_soln(a, b, c):
    #returns the smallest positive value that solves the quadratic equation

    t1 = (-b + math.sqrt(b**2 - 4*a*c))/(2*a)
    t2 = (-b - math.sqrt(b**2 - 4*a*c))/(2*a)

    if (abs(t1) < 10**-8):
        t1 = 0
    if (abs(t2) < 10**-8):
        t2 = 0

    if (abs(t1)-abs(t2)) > 10**83:

        if t1 < t2:
            t1 = 0
        elif t2 < t1:
            t2 = 0

    if t1 > 0 and t2 > 0:
        return min(t1,t2)
    elif t1 > 0 and t2 <= 0:
        return t1
    elif t1 <= 0 and t2 > 0:
        return t2
    else:
        raise Exception("Error: no positive collision times.")


def compute_collision_time(s1, s2):
    #only consider sphere pairs that are approaching each other, aka if <(r1-r2),(v1-v2)> less than 0
    pos_diff_x = s2.pos_x-s1.pos_x
    pos_diff_y = s2.pos_y-s1.pos_y
    pos_diff_z = s2.pos_z-s1.pos_z

    vel_diff_x = s2.vel_x-s1.vel_x
    vel_diff_y = s2.vel_y-s1.vel_y
    vel_diff_z = s2.vel_z-s1.vel_z

    #calculate quadratic equation coefficients based on current sphere pair
    ## a = <v2-v1,v2-v1>
    a = (vel_diff_x**2 + vel_diff_y**2 + vel_diff_z**2)
    ## b = <2(p2-p1),v2-v1>
    b = 2*(pos_diff_x*vel_diff_x + pos_diff_y*vel_diff_y + pos_diff_z*vel_diff_z)
    ## c = <p2-p1.p2-p1> - (radius1 + radius2)**2
    c = (pos_diff_x**2 + pos_diff_y**2 + pos_diff_z**2) - (s1.radius + s2.radius)**2

    if (pos_diff_x*vel_diff_x + pos_diff_y*vel_diff_y + pos_diff_z*vel_diff_z) < 0 and (b**2-4*a*c) >= 0:
        return minimum_pos_quadratic_soln(a,b,c)
    else:
        return -2


def compute_reflection_time(sphere, universe_radius):
    #check when each sphere would hit the wall if moving and not blocked

    ##checking to see if sphere is moving
    if not(sphere.vel_x == 0 and sphere.vel_y == 0 and sphere.vel_z == 0):

        a = sphere.vel_x**2 + sphere.vel_y**2 + sphere.vel_z**2

        b = (2*sphere.pos_x)*sphere.vel_x + (2*sphere.pos_y)*sphere.vel_y + (2*sphere.pos_z)*sphere.vel_z

        c = sphere.pos_x**2 + sphere.pos_y**2 + sphere.pos_z**2 - (universe_radius - sphere.radius)**2

        return minimum_pos_quadratic_soln(a,b,c)

    return -2

def compute_positions(sphere_list, event_time):
    for sphere in sphere_list:
        sphere.pos_x += sphere.vel_x*event_time
        sphere.pos_y += sphere.vel_y*event_time
        sphere.pos_z += sphere.vel_z*event_time

# get collsion results
def compute_collision(sphere_a, sphere_b):
	# compile variables
    m1 = sphere_a.mass
    m2 = sphere_b.mass

    r1 = sphere_a.radius
    r2 = sphere_b.radius

    pos1_x = sphere_a.pos_x
    pos1_y = sphere_a.pos_y
    pos1_z = sphere_a.pos_z
    pos2_x = sphere_b.pos_x
    pos2_y = sphere_b.pos_y
    pos2_z = sphere_b.pos_z

    v1_x = sphere_a.vel_x
    v1_y = sphere_a.vel_y
    v1_z = sphere_a.vel_z
    v2_x = sphere_b.vel_x
    v2_y = sphere_b.vel_y
    v2_z = sphere_b.vel_z

    # compute inbetween values
    dot_product_1_top = (v1_x-v2_x)*(pos1_x-pos2_x) + (v1_y-v2_y)*(pos1_y-pos2_y) + (v1_z-v2_z)*(pos1_z-pos2_z)
    dot_product_1_bottom = (pos1_x-pos2_x)**2 + (pos1_y-pos2_y)**2 + (pos1_z-pos2_z)**2

    dot_product_2_top = (v2_x-v1_x)*(pos2_x-pos1_x) + (v2_y-v1_y)*(pos2_y-pos1_y) + (v2_z-v1_z)*(pos2_z-pos1_z)
    dot_product_2_bottom = (pos2_x-pos1_x)**2 + (pos2_y-pos1_y)**2 + (pos2_z-pos1_z)**2

    # compute new velocity values for each colliding sphere
    v1_x_prime = v1_x - (2*m2*dot_product_1_top*(pos1_x-pos2_x)) / ((m1+m2)*dot_product_1_bottom)
    v1_y_prime = v1_y - (2*m2*dot_product_1_top*(pos1_y-pos2_y)) / ((m1+m2)*dot_product_1_bottom)
    v1_z_prime = v1_z - (2*m2*dot_product_1_top*(pos1_z-pos2_z)) / ((m1+m2)*dot_product_1_bottom)

    v2_x_prime = v2_x - (2*m1*dot_product_2_top*(pos2_x-pos1_x)) / ((m1+m2)*dot_product_2_bottom)
    v2_y_prime = v2_y - (2*m1*dot_product_2_top*(pos2_y-pos1_y)) / ((m1+m2)*dot_product_2_bottom)
    v2_z_prime = v2_z - (2*m1*dot_product_2_top*(pos2_z-pos1_z)) / ((m1+m2)*dot_product_2_bottom)

    sphere_a.vel_x = v1_x_prime
    sphere_a.vel_y = v1_y_prime
    sphere_a.vel_z = v1_z_prime

    sphere_b.vel_x = v2_x_prime
    sphere_b.vel_y = v2_y_prime
    sphere_b.vel_z = v2_z_prime


# get reflection results
def compute_reflection(sphere):

    pos_x = sphere.pos_x
    pos_y = sphere.pos_y
    pos_z = sphere.pos_z

    v_x = sphere.vel_x
    v_y = sphere.vel_y
    v_z = sphere.vel_z

    dot_product_top = (v_x)*(pos_x) + (v_y)*(pos_y) + (v_z)*(pos_z)
    dot_product_bottom = (pos_x)**2 + (pos_y)**2 + (pos_z)**2

    sphere.vel_x -= 2*(dot_product_top/dot_product_bottom)*pos_x
    sphere.vel_y -= 2*(dot_product_top/dot_product_bottom)*pos_y
    sphere.vel_z -= 2*(dot_product_top/dot_product_bottom)*pos_z

def compute_energy(sphere_list):
    sum_energy = 0
    for sphere in sphere_list:
        sum_energy += (0.5 * sphere.mass * (sphere.vel_x**2+sphere.vel_y**2+sphere.vel_z**2))
    return sum_energy


def compute_momentum(sphere_list):

    momentum_x = 0
    momentum_y = 0
    momentum_z = 0
    for sphere in sphere_list:
        momentum_x += sphere.mass*sphere.vel_x
        momentum_y += sphere.mass*sphere.vel_y
        momentum_z += sphere.mass*sphere.vel_z

    moment = [momentum_x,momentum_y,momentum_z]
    return moment

# Our main
def run_sim(universe_radius, duration):
    for string in sphere_strings:
        sphere_list.append(Sphere(string))

    #print out all initial conditions
    print("\nHere are the initial conditions.")
    print(f"universe radius {universe_radius:g}")
    print(f"end simulation {duration:g}")

    for sphere in sphere_list:
        print(f"{sphere.name} m={sphere.mass:g} R={sphere.radius:g} p=({sphere.pos_x:g},{sphere.pos_y:g},{sphere.pos_z:g}) v=({sphere.vel_x:g},{sphere.vel_y:g},{sphere.vel_z:g})")

    #calculate initial energy and momentum
    energy = compute_energy(sphere_list)
    momentum = compute_momentum(sphere_list)

    print(f"energy: {energy:g}")
    print(f"momentum: ({momentum[0]:g},{momentum[1]:g},{momentum[2]:g})")
    # begin simulation
    time_elapsed = 0

    #end simulation once the nearest event would occur beyond existence of universe
    while time_elapsed < duration:

        simul_collisions=[]
        checked =0
        nearest_event_time = -1
        ## determine the time of the first event by evaluating all possible collision times from starting
        # positions, and choosing event that happens soonest

        next_colliding_pair_indices = []
        # check sphere-to-sphere collision times
        for i in range(len(sphere_list)-1):
            s1 = sphere_list[i]
            j = i+1

            while j < len(sphere_list):
                s2 = sphere_list[j]

                current_event_time = compute_collision_time(s1,s2)
                # print("current event time: " + str(current_event_time))
                # print("nearest event time: " + str(nearest_event_time))

                #print("collision current event time:",current_event_time) DEBUGGER
                # compare current event time to value of nearest_event_time to see if this event would happen sooner, and update nearest_event_time if so
                if (current_event_time < nearest_event_time or nearest_event_time == -1) and current_event_time != -2 and current_event_time != 0:
                    nearest_event_time = current_event_time
                    next_event_type = "colliding"
                    next_colliding_pair_indices = [i, j]
                    simul_collisions = []
                    simul_collisions.append([s1,s2])
                    #print("initial len"+str(len(simul_collisions)))
                    checked = 1

                if(current_event_time == nearest_event_time and nearest_event_time > 0) and current_event_time != -2 and current_event_time != 0:
                    if([s1,s2] not in simul_collisions):
                        #print(len(simul_collisions))
                        simul_collisions.append([s1,s2])

                # increment j to look at next possible sphere pair
                j += 1


        # check sphere-to-wall collision times
        for sphere in sphere_list:
            # check when each sphere would next hit the wall if moving and not blocked

            current_event_time = compute_reflection_time(sphere,universe_radius)

            # compare current event time to value of nearest_event_time to see if this event would happen sooner, and update nearest_event_time if so
            if (current_event_time < nearest_event_time or nearest_event_time == -1) and current_event_time != -2 and current_event_time != 0:
                nearest_event_time = current_event_time
                next_reflecting_sphere = sphere
                next_event_type = "reflecting"
                checked = 1

        ## DEBUG:
        # print("Collision pairs:")
        # for pair in simul_collisions:
        #     print(f"{pair[0].name}, {pair[1].name}")

        #make adjustments to sphere velocities based on next occurring event (use nearest_event_time and next_event_type), and add event to event list using event class
        #print("\nnearest event time:",nearest_event_time) #DEBUGGER
        #print("nearest event type:",next_event_type) #DEBUGGER
        time_elapsed += nearest_event_time
        # print("time_elapsed:",time_elapsed)

        if time_elapsed < duration:
            #compute new positions for each of the spheres based on event time
            compute_positions(sphere_list, nearest_event_time)
            multi_collisions = 0
            #print(len(simul_collisions))
            #print(simul_collisions)
            if len(simul_collisions) > 2:
                multi_collisions = 1

                primary_sphere_touch_count = 0
                previous_primary_sphere = "NA"

                for pairs in simul_collisions:

                    # Compute new velocities, system energy, and system momentum
                    colliding_sphere1 = pairs[0]
                    colliding_sphere2 = pairs[1]

                    current_primary_sphere = colliding_sphere1.name

                    if current_primary_sphere == previous_primary_sphere:
                        primary_sphere_touch_count += 1
                    else:
                        primary_sphere_touch_count = 0

                    if primary_sphere_touch_count < 2:
                        compute_collision(colliding_sphere1, colliding_sphere2)

                        updated_energy = compute_energy(sphere_list)
                        updated_momentum = compute_momentum(sphere_list)

                        # create an event and append to event list
                        current_event = Event(time_elapsed, next_event_type, colliding_sphere1, colliding_sphere2, sphere_list, updated_energy, updated_momentum)
                        event_list.append(current_event)

                    previous_primary_sphere = current_primary_sphere

                simul_collisions = []

            if next_event_type == "reflecting" and multi_collisions == 0:
                # Compute new velocities, system energy, and system momentum
                compute_reflection(next_reflecting_sphere)
                updated_energy = compute_energy(sphere_list)
                updated_momentum = compute_momentum(sphere_list)

                #create an event and append to event list
                current_event = Event(time_elapsed, next_event_type, next_reflecting_sphere, None, sphere_list, energy, updated_momentum)
                event_list.append(current_event)

            elif next_event_type == "colliding" and multi_collisions == 0:
                # Compute new velocities, system energy, and system momentum
                colliding_sphere1 = sphere_list[next_colliding_pair_indices[0]]
                colliding_sphere2 = sphere_list[next_colliding_pair_indices[1]]

                compute_collision(colliding_sphere1, colliding_sphere2)

                updated_energy = compute_energy(sphere_list)
                updated_momentum = compute_momentum(sphere_list)

                # create an event and append to event list
                current_event = Event(time_elapsed, next_event_type, colliding_sphere1, colliding_sphere2, sphere_list, updated_energy, updated_momentum)
                event_list.append(current_event)

            elif multi_collisions==0 and (next_event_type != "reflecting" or next_event_type != "colliding"):
                raise Exception("Error: no event type specified")

            simul_collisions = []

            # #DEBUGGER: ERASE AFTER 100% SUCCESSFUL TESTING
            # for sphere in sphere_list:
            #     print(f"{sphere.name} p=({sphere.pos_x:g},{sphere.pos_y:g},{sphere.pos_z:g}) v=({sphere.vel_x:g},{sphere.vel_y:g},{sphere.vel_z:g})")
            #
            # print("\n")

    print("\nHere are the events.\n")

    for event in event_list:
        report_event(event)

    #terminate simulation after successful error report
    sys.exit(0)


def report_event(event):
    print(f"time of event: {event.timestamp:g}")

    if event.type == "colliding":
        print(event.type, event.sphere_names_involved[0], event.sphere_names_involved[1])
    else:
        print(event.type, event.sphere_names_involved[0])

    for log in event.sphere_status_logs:
        print(log)

    print(f"energy: {event.energy:g}")
    print(f"momentum: ({event.momentum[0]:g},{event.momentum[1]:g},{event.momentum[2]:g})\n")

######################################################################################################
######################################################################################################
######################################################################################################

# Script starts here
sphere_strings = []
sphere_list = []
event_list = []
simul_collisions = []
simul_reflections =[]

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
