import sys
import numpy as np
from scipy.integrate import odeint
import argparse
# The gravitational acceleration (m.s-2).
g = 9.81


def deriv(y, t, L1, L2, m1, m2):
    """Return the first derivatives of y = theta1, z1, theta2, z2."""
    theta1, z1, theta2, z2 = y

    c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)

    theta1dot = z1
    z1dot = (m2*g*np.sin(theta2)*c - m2*s*(L1*z1**2*c + L2*z2**2) -
             (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
    theta2dot = z2
    z2dot = ((m1+m2)*(L1*z1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) +
             m2*L2*z2**2*s*c) / L2 / (m1 + m2*s**2)
    return theta1dot, z1dot, theta2dot, z2dot

def solve(L1, L2, m1, m2, tmax, dt, y0):
    t = np.arange(0, tmax+dt, dt)

    # Do the numerical integration of the equations of motion
    y = odeint(deriv, y0, t, args=(L1, L2, m1, m2))
    theta1, theta2 = y[:,0], y[:,2]

    # Convert to Cartesian coordinates of the two bob positions.
    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)
    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)

    return theta1, theta2, x1, y1, x2, y2

def simulate_pendulum(theta_resolution,tmax,dt,resoutls_file):
    # Pendulum rod lengths (m), bob masses (kg).
    L1, L2 = 1.0, 1.0
    m1, m2 = 1.0, 1.0

    # Maximum time, time point spacings (all in s).
    # tmax, dt = 30.0, 0.01

    # Means to run to csv
    import csv
    with open(resoutls_file+'.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(['theta1_init', 'theta2_init',
                                     'theta1',
                                     'theta2',
                                     'x1',
                                     'y1',
                                     'x2',
                                     'y2'])

        for theta1_init in np.linspace(0, 2*np.pi, theta_resolution):
            for theta2_init in np.linspace(0, 2*np.pi, theta_resolution):
                # Initial conditions: theta1, dtheta1/dt, theta2, dtheta2/dt.
                y0 = np.array([
                    theta1_init,
                    0.0,
                    theta2_init,
                    0.0
                ])

                theta1, theta2, x1, y1, x2, y2 = solve(L1, L2, m1, m2, tmax, dt, y0)
                spamwriter.writerow([theta1_init, theta2_init,
                                     theta1[-1],
                                     theta2[-1],
                                     x1[-1],
                                     y1[-1],
                                     x2[-1],
                                     y2[-1]])
                # print theta1_init, theta2_init, theta1[-1], theta2[-1]

def main():
    parser = argparse.ArgumentParser(
    )
    parser.add_argument(
        '-r',
        '--resolution',
        dest='resolution',
        metavar='resolution',
        default=5,
        help="Define lvl of precision"
    )
    parser.add_argument(
        '-dt',
        '--delta-time',
        dest='dt',
        metavar='dt',
        default=0.01,
        help="Define dt"
    )
    parser.add_argument(
        '-tm',
        '--time-max',
        dest='time_max',
        metavar='time_max',
        default=30,
        help="Define lvl of precision"
    )
    parser.add_argument(
        '-f',
        '--results-file',
        dest='results_file',
        metavar='results_file',
        default='results',
        help="Define name of results file"
    )
    args = parser.parse_args()

    resolution = float(args.resolution)
    resoutls_file = args.results_file
    dt = float(args.dt)
    tmax = float(args.time_max)
    print "resolution: ", resolution
    print "dt: ", dt
    print "tmax: ", tmax
    print "result file: ", resoutls_file

    simulate_pendulum(resolution,tmax,dt,resoutls_file)

if __name__ == "__main__":
    main()

