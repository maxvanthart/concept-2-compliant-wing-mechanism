"""Reference geometry for a user-specified planar PRBM."""
import numpy as np


# One angle per rigid body. L4 and L9 set the axes of the two triangles.
FLEXURE_PAIRS = {
    'B_BQ': (0, 5), 'B_BD': (0, 10), 'P': (1, 4), 'Q': (5, 4),
    'D': (10, 9), 'T1': (4, 11), 'T2': (9, 11),
}


def relative_flexure_angles(positions, a, b):
    """Continuous relative angle of body b to body a along the pose sequence."""
    return np.unwrap([pos[f'theta{b}'] - pos[f'theta{a}'] for pos in positions])


def reference_geometry(virtual_springs, rigid_links, fixed_points=None):
    """Evaluate virtual pivots and pivot-to-pivot lengths in one reference pose.

    Each spring supplies two CAD points on its reference flexure axis in metres
    and a model-supplied pivot fraction measured from the first point.
    Each rigid link supplies the names of its two pivot or fixed reference points.
    """
    def point(coordinates):
        xy = np.asarray(coordinates, dtype=float)
        if xy.shape != (2,) or not np.all(np.isfinite(xy)):
            raise ValueError('Reference coordinates must be finite planar points in metres')
        return xy

    points = {name: point(xy) for name, xy in (fixed_points or {}).items()}
    springs = {}
    for name, spec in virtual_springs.items():
        start = point(spec['start_xy_m'])
        end = point(spec['end_xy_m'])
        fraction = float(spec['pivot_fraction'])
        if name in points or not np.any(end != start) or not 0 <= fraction <= 1:
            raise ValueError(f'{name}: expected a unique pivot, distinct endpoints and fraction in [0, 1]')
        points[name] = start + fraction * (end - start)
        springs[name] = points[name].tolist()

    lengths = {}
    for link, (a, b) in rigid_links.items():
        if a not in points or b not in points:
            raise ValueError(f'{link}: unknown pivot or fixed reference point')
        length = float(np.linalg.norm(points[b] - points[a]))
        if length <= 0:
            raise ValueError(f'{link}: effective rigid length must be positive')
        lengths[link] = length
    return {'virtual_spring_xy_m': springs, 'effective_link_length_m': lengths}
