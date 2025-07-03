import math
import numpy as np

def tangential_velocity(n, d):
    """ Convert rpm and diameter to tangential velocity

    :param n (float): rpm [1/min]
    :param d (float): diameter [m]
    :return u (float): tangential velocity [m/s]
    """
    u = np.pi*d*n/60
    return u

def interpolate_curve(points, num_points=200):
    """
    Given a list of points (tuples, e.g., [(x1, y1), (x2, y2), ...]) along a curve,
    interpolate the curve to produce a new list of points that are equally spaced 
    in arc-length.
    
    Parameters:
        points (list of tuples): Original points along the curve.
        num_points (int): Desired number of points in the output.
        
    Returns:
        new_points (list of tuples): Interpolated points with equal arc-length spacing.
    """
    if len(points) < 2:
        return points

    # Compute cumulative distances along the curve.
    distances = [0.0]
    for i in range(1, len(points)):
        dx = points[i][0] - points[i-1][0]
        dy = points[i][1] - points[i-1][1]
        d = math.hypot(dx, dy)
        distances.append(distances[-1] + d)
    total_length = distances[-1]
    
    # Create equally spaced distances along the total length.
    new_distances = np.linspace(0, total_length, num_points)
    
    new_points = []
    # For each target distance, find the segment in which it falls.
    for nd in new_distances:
        # Find index i such that distances[i] <= nd <= distances[i+1]
        for i in range(len(distances) - 1):
            if distances[i] <= nd <= distances[i+1]:
                # Linear interpolation factor
                segment_length = distances[i+1] - distances[i]
                t = (nd - distances[i]) / segment_length if segment_length != 0 else 0
                x = points[i][0] + t * (points[i+1][0] - points[i][0])
                y = points[i][1] + t * (points[i+1][1] - points[i][1])
                new_points.append((x, y))
                break
    return new_points