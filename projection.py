# モルワイデ図法の座標変換については別紙。
# https://colab.research.google.com/drive/1XAcSpD2fKkiu3o69CJ93okSR_-Rput3A?usp=sharing

# 経度緯度と、極座標では緯度(縦方法)の表しかたが違うことに注意。経度緯度では緯度0は赤道だが、極座標では$\theta=0$は北極を表す。
import numpy as np

def equirectangular_grid(width=1000):
    """
    The function `equirectangular_grid` generates a grid of points in an equirectangular projection with
    a specified width.
    
    Args:
      width: The width parameter determines the number of points along the x-axis of the grid. It is
    used to create a grid with a specified width and height. Defaults to 1000
    
    Returns:
      a meshgrid of x and y values.
    """
    height = width // 2

    x = np.linspace(-np.pi, +np.pi, width)
    y = np.linspace(-np.pi/2, +np.pi/2, height)

    return np.meshgrid(x, y)


def polar(x,y,z)->tuple:
    """
    The function `polar` takes three input parameters `x`, `y`, and `z` and returns a tuple containing
    the polar angles of the vector (x, y, z).
    
    Args:
      x: The x-coordinate of the point in Cartesian coordinates.
      y: The y-coordinate in a Cartesian coordinate system.
      z: The parameter `z` represents the height or vertical component in a three-dimensional coordinate
    system.
    
    Returns:
      a tuple containing two values. The first value is the arctangent of y/x, and the second value is
    the arctangent of z divided by the square root of (x^2 + y^2).
    """
    return np.arctan2(y,x), np.arctan2(z,(x**2+y**2)**0.5)


def mollweide_grid(width=1000):
    """
    The function `mollweide_grid` returns the longitude `alpha` and latitude `delta` of a point on a
    Mollweide projection given its pixel coordinates `(x, y)`.
    
    Args:
      width: the width of the grid in pixels. It determines the resolution of
    the Mollweide projection. Defaults to 1000
    
    Returns:
      two arrays: `alpha` and `delta`. Both `alpha` and `delta` are 1000x500 two-dimensional arrays. `alpha` and `beta` represents the longitude and latitude (in radians) of each pixel (x,
    y) in the Mollweide projection. `alpha` is None when it is out of the boundary.
    """
    X, Y = equirectangular_grid(width)

    theta = np.arcsin(2*Y/np.pi)
    alpha = X / ((1-(2*Y/np.pi)**2)**0.5)
    alpha[alpha<-np.pi] = np.nan
    alpha[alpha>np.pi] = np.nan
    delta = np.arcsin((2*theta+np.sin(2*theta)) / np.pi)

    return alpha, delta