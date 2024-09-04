import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import csv

G = 6.67 * 10 ** (-11)
C = 0.00596
AU = 1.496 * (10 ** (11))
time_step = 1000000
print("The following simulation assumes point-like masses")
print("Distances in AU and mass in solar-mass")

filename = "parameters.csv"


class simulation_state:
    def __init__(self, filename):
        (
            self.body_id,
            self.mass,
            self.x_pos,
            self.y_pos,
            self.z_pos,
            self.vX,
            self.vY,
            self.vZ,
            self.aX,
            self.aY,
            self.aZ,
        ) = read_parameters(filename)
        self.time_step = 1000000
        self.x_positions = np.array(self.x_pos)
        self.y_positions = np.array(self.y_pos)
        self.z_positions = np.array(self.z_pos)


def read_parameters(filename):
    body_ids = []
    masses = []
    xs = []
    ys = []
    zs = []
    vXs = []
    vYs = []
    vZs = []
    aXs = []
    aYs = []
    aZs = []

    try:
        with open(filename, mode="r") as file:
            csv_reader = csv.DictReader(file)

            for row in csv_reader:
                body_ids.append(int(row["body_id"]))
                masses.append(float(row["mass"]))
                xs.append(float(row["x"]))
                ys.append(float(row["y"]))
                zs.append(float(row["z"]))
                vXs.append(float(row["vx"]))
                vYs.append(float(row["vy"]))
                vZs.append(float(row["vz"]))
                aXs.append(float(row["ax"]))
                aYs.append(float(row["ay"]))
                aZs.append(float(row["az"]))

    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")
        return None
    except KeyError as e:
        print(f"Error: Missing column {e} in the CSV file.")
        return None
    except ValueError as e:
        print(f"Error: Invalid value found in the CSV file. {e}")
        return None

    return (
        np.array(body_ids),
        np.array(masses),
        np.array(xs),
        np.array(ys),
        np.array(zs),
        np.array(vXs),
        np.array(vYs),
        np.array(vZs),
        np.array(aXs),
        np.array(aYs),
        np.array(aZs),
    )


state = simulation_state(filename)
n = len(state.body_id)


# Get acceleration
def findaccln(x_pos, y_pos, z_pos, aX, aY, aZ, mass):
    for i in range(n):
        for j in range(n):
            if i > j:
                dist = (
                    ((x_pos[i] - x_pos[j]) ** 2)
                    + ((y_pos[i] - y_pos[j]) ** 2)
                    + ((z_pos[i] - z_pos[j]) ** 2)
                ) ** (3 / 2)
                aX[i] = -(C * mass[j] * (x_pos[i] - x_pos[j])) / dist
                aY[i] = -(C * mass[j] * (y_pos[i] - y_pos[j])) / dist
                aZ[i] = -(C * mass[j] * (z_pos[i] - z_pos[j])) / dist
                aX[j] = -(mass[i] / mass[j]) * (aX[i])
                aY[j] = -(mass[i] / mass[j]) * (aY[i])
                aZ[j] = -(mass[i] / mass[j]) * (aZ[i])
    return aX, aY, aZ


# Update
def updateproperties(x_pos, y_pos, z_pos, vX, vY, vZ, aX, aY, aZ, mass, timestep):
    aX, aY, aZ = findaccln(x_pos, y_pos, z_pos, aX, aY, aZ, mass)
    vX = vX + aX * timestep
    vY = vY + aY * timestep
    vZ = vZ + aZ * timestep
    x_pos = x_pos + vX * timestep / (AU)
    y_pos = y_pos + vY * timestep / (AU)
    z_pos = z_pos + vZ * timestep / (AU)
    return x_pos, y_pos, z_pos, vX, vY, vZ, aX, aY, aZ


Kinetic_Energy = [[0] for i in range(n)]
Potential_Energy = [[0.0 for j in range(n)] for k in range(n)]


fig = plt.figure()
ax = fig.add_subplot(projection="3d")


ax.set_xlim(-20, 20)
ax.set_ylim(-20, 20)
ax.set_zlim(-20, 20)

scatter = ax.scatter([], [], [])


def update(frame):
    (
        state.x_pos,
        state.y_pos,
        state.z_pos,
        state.vX,
        state.vY,
        state.vZ,
        state.aX,
        state.aY,
        state.aZ,
    ) = updateproperties(
        state.x_pos,
        state.y_pos,
        state.z_pos,
        state.vX,
        state.vY,
        state.vZ,
        state.aX,
        state.aY,
        state.aZ,
        state.mass,
        state.time_step,
    )

    state.x_positions = np.append(state.x_positions, state.x_pos)
    state.y_positions = np.append(state.y_positions, state.y_pos)
    state.z_positions = np.append(state.z_positions, state.z_pos)

    scatter._offsets3d = (state.x_pos, state.y_pos, state.z_pos)
    return (scatter,)


# Create animation
ani = animation.FuncAnimation(
    fig, update, frames=1000, interval=20, blit=False, repeat=True
)


# Set labels and title
ax.set_xlabel("X axis")
ax.set_ylabel("Y axis")
ax.set_zlabel("Z axis")
ax.set_title("3D Visualization of Gravitational Bodies")

# Show the plot
plt.show()

