import numpy as np
from os.path import exists

G = 6.67430 * 10 ** (-11)

number = int(input("Select how many bodies you want: "))
time_step = float(
    input("Input time step. A smaller time step gives more accurate results: ")
)
lengthofsimul = float(
    input("Input the amount of time you want to run the simulation for: ")
)
print("The following simulation assumes point-like masses")
print("Distances in AU and mass in solar-mass")

# Initialize Measurable quantities
abs_pos = [[0.0, 0.0, 0.0] for i in range(number)]
rel_pos = [[[0.0, 0.0, 0.0] for j in range(number)] for i in range(number)]
velocity = [[0.0, 0.0, 0.0] for i in range(number)]
mass = [[0.0] for i in range(number)]
accln = [[0.0, 0.0, 0.0] for i in range(number)]

# Take Values from file
if exists("nbodyparameters.csv"):
    f = open("nbodyparameters", "r")
    pass


# Accept initial parameters
def AcceptParameters(position, velocity, mass):
    for i in range(number):
        # Position
        position[i][0] = float(
            input("Please input Co-ordinate {} of particle {}: ".format(1, i + 1))
        )
        position[i][1] = float(
            input("Please input Co-ordinate {} of particle {}: ".format(2, i + 1))
        )
        position[i][2] = float(
            input("Please input Co-ordinate {} of particle {}: ".format(3, i + 1))
        )
        # Velocity
        velocity[i][0] = float(
            input(
                "Please input Velocity Component {} of particle {}: ".format(1, i + 1)
            )
        )
        velocity[i][1] = float(
            input(
                "Please input Velocity Component {} of particle {}: ".format(2, i + 1)
            )
        )
        velocity[i][2] = float(
            input(
                "Please input Velocity Component {} of particle {}: ".format(3, i + 1)
            )
        )
        # Mass
        mass[i] = float(input("Please input Mass of particle {}: ".format(i + 1)))
    return position, velocity, mass


# Get acceleration
def findaccln(position, ds, accln):
    for i in range(number):
        accln[i][0] = 0
        accln[i][1] = 0
        accln[i][2] = 0
        for j in range(number):
            if i > j:
                # Find Relative Position
                ds[i][j][0] = position[i][0] - position[j][0]
                ds[i][j][1] = position[i][1] - position[j][1]
                ds[i][j][2] = position[i][2] - position[j][2]
                dist = (
                    ((ds[i][j][0]) ** 2) + ((ds[i][j][1]) ** 2) + ((ds[i][j][2]) ** 2)
                ) ** (3 / 2)
                accx = -(0.00596 * mass[j] * ds[i][j][0]) / dist
                accy = -(0.00596 * mass[j] * ds[i][j][1]) / dist
                accz = -(0.00596 * mass[j] * ds[i][j][2]) / dist
                accln[i][0] = accln[i][0] + (accx)
                accln[i][1] = accln[i][1] + (accy)
                accln[i][2] = accln[i][2] + (accz)
                accln[j][0] = accln[j][0] - (mass[i] / mass[j]) * (accx)
                accln[j][1] = accln[j][1] - (mass[i] / mass[j]) * (accy)
                accln[j][2] = accln[j][2] - (mass[i] / mass[j]) * (accz)
    return ds, accln


abs_pos, velocity, mass = AcceptParameters(abs_pos, velocity, mass)
rel_pos, accln = findaccln(abs_pos, rel_pos, accln)


# Update
def updateproperties(position, velocity, accln, timestep, ds):
    ds, accln = findaccln(abs_pos, rel_pos, accln)
    # print(accln)
    for i in range(number):
        velocity[i][0] = velocity[i][0] + accln[i][0] * timestep
        velocity[i][1] = velocity[i][1] + accln[i][1] * timestep
        velocity[i][2] = velocity[i][2] + accln[i][2] * timestep
        position[i][0] = position[i][0] + velocity[i][0] * timestep / (
            1.496 * (10 ** (11))
        )
        position[i][1] = position[i][1] + velocity[i][1] * timestep / (
            1.496 * (10 ** (11))
        )
        position[i][2] = position[i][2] + velocity[i][2] * timestep / (
            1.496 * (10 ** (11))
        )
    return position, velocity, accln, ds


Kinetic_Energy = [[0] for i in range(number)]
Potential_Energy = [[0.0 for j in range(number)] for k in range(number)]


# Calculate Energies
def energy(position, ds, velocity, mass, G):
    for i in range(number):
        # Kinetic energy here
        Kinetic_Energy[i] = (
            (
                0.5
                * mass[i]
                * (velocity[i][0] ** 2 + velocity[i][1] ** 2 + velocity[i][2] ** 2)
            )
            * 1.989
            * (10**30)
        )
        for j in range(number):
            if i > j:
                # Potential energy here

                ds[i][j][0] = position[i][0] - position[j][0]
                ds[i][j][1] = position[i][1] - position[j][1]
                ds[i][j][2] = position[i][2] - position[j][2]
                dist = (
                    ((ds[i][j][0]) ** 2) + ((ds[i][j][1]) ** 2) + ((ds[i][j][2]) ** 2)
                ) ** (1 / 2)
                Potential_Energy[i][j] = -(
                    G * mass[i] * mass[j] * ((1.989 * (10**30)) ** 2) / dist
                )
    total_energy = sum(Kinetic_Energy) + sum(sum(Potential_Energy, []))
    return total_energy


# Display properties
for k in range(round(lengthofsimul / time_step)):
    abs_pos, velocity, accln, rel_pos = updateproperties(
        abs_pos, velocity, accln, time_step, rel_pos
    )
    print(abs_pos)
    # # Use To verify the code
    # total_energy = energy(abs_pos, rel_pos, velocity, mass, G)
    # print(total_energy)
