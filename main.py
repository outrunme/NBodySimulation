import numpy as np

number = int(input("Select how many bodies you want: "))
time_step = float(
    input("Input time step. A smaller time step gives more accurate results: ")
)
lengthofsimul = float(
    input("Input the amount of time you want to run the simulation for: ")
)
print("The following simulation assumes point-like masses")
print("Distances in AU and mass in solar-mass")
abs_pos = [[0.0 for j in range(2)] for i in range(number)]
rel_pos = [[[0.0 for k in range(2)] for j in range(number)] for i in range(number)]
velocity = [[0.0 for j in range(2)] for i in range(number)]
acceleration = [[0.0 for j in range(2)] for i in range(number)]
mass = [0] * number


def AcceptPosition(position):
    for i in range(number):
        for j in range(2):
            position[i][j] = float(
                input(
                    "Please input Co-ordinate {} of particle {}: ".format(j + 1, i + 1)
                )
            )
    return position


def AcceptVelocity(velocity):
    for i in range(number):
        for j in range(2):
            velocity[i][j] = float(
                input(
                    "Please input Velocity Component {} of particle {}: ".format(
                        j + 1, i + 1
                    )
                )
            )
    return velocity


def AcceptMass(mass):
    for i in range(number):
        print("PLease input Mass of particle", (i + 1))
        mass[i] = float(input())
    return mass


def definerelativeposition(relativepostion, position):
    for i in range(number):
        for j in range(number):
            for k in range(2):
                relativepostion[i][j][k] = position[i][k] - position[j][k]
    return relativepostion


def findacceleration(relativeposition, acceleration):
    for i in range(number):
        for j in range(number):
            for k in range(2):
                if not i == j:
                    acceleration[i][k] = acceleration[i][k] + ((0.005929 * mass[j])) * (
                        -relativeposition[i][j][k]
                    ) / (
                        (relativeposition[i][j][0]) ** 2
                        + ((relativeposition[i][j][1]) ** 2)
                    ) ** (
                        3 / 2
                    )
    return acceleration


abs_pos = AcceptPosition(abs_pos)
velocity = AcceptVelocity(velocity)
mass = AcceptMass(mass)
rel_pos = definerelativeposition(rel_pos, abs_pos)
acceleration = findacceleration(rel_pos, acceleration)


# update
def updateproperties(position, velocity, acceleration, timestep, relativeposition):
    for i in range(number):
        for j in range(2):
            position[i][j] = position[i][j] + velocity[i][j] * timestep
            velocity[i][j] = velocity[i][j] + acceleration[i][j] * timestep
            acceleration = findacceleration(relativeposition, acceleration)
    return position, velocity, acceleration


# abs_pos, velocity = updateproperties(
#     abs_pos, velocity, acceleration, time_step, rel_pos
# )
for k in range(round(lengthofsimul / time_step)):
    abs_pos, velocity, acceleration = updateproperties(
        abs_pos, velocity, acceleration, time_step, rel_pos
    )
    print(abs_pos, velocity, acceleration)

print(abs_pos)
print(velocity)
print(acceleration)
