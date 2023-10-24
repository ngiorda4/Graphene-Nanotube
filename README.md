# Graphene-Nanotube
Build a graphene sheet and roll it up in order to construct a carbon nanotube!

## Description
In this project, a graphene sheet was constructed and rolled up to form a carbon nanotube. For the graphene sheet, orientation was determined with the origin set at the atom located at the leftmost corner of the sheet. For the nanotube, the z-axis passed through the center of the tube, and the x-axis through the atom originally at the origin. All lengths were measured in units of the carbon-carbon bond length, meaning bonded atoms were a distance of 1.0 apart.

To visually represent the graphene sheet and its corresponding carbon nanotube, four functions were coded within a Python script named `nano_plot.py`. These functions visually represented atoms as circles (for a graphene sheet) and depicted carbon-carbon bonds as lines.

## Functions in `nano_plot.py`
1. **find_connectivity**: This function calculates the Euclidean distance squared between two points to determine their connectivity. The distance is compared against an upper and lower threshold, and if it falls outside this range, its value is set to zero.

**Parameters:**
- **points**: A 2D numpy array containing x and y coordinates of points in 2D space. The dimensions of this numpy array are (number of points × 2).
- **low_th**: A scalar float defining the lower threshold for the Euclidean distance squared. Default value: 0.95.
- **up_th**: A scalar float defining the upper threshold for the Euclidean distance squared. Default value: 1.05.

**Returns:**  
A 2D numpy array where elements have a value of 0 if the Euclidean distance squared between the i-th and j-th points is ≤ 0.95 or ≥ 1.05. Otherwise, the value is computed as:  
For points \( p_1 = (x_1, y_1) \) and \( p_2 = (x_2, y_2) \), the Euclidean distance squared is \( (x_2 - x_1)^2 + (y_2 - y_1)^2 \).

2. **get_connectivity_ind**: The second function, get connectivity ind, searches the distance numpy array and returns the indexes of the elements that contain a Euclidean distance squared value greater than 0. In this function, we will call the find connectivity function and then search for the indexes of the relevant connected points in the distance array.

**Parameters:**
- **points**: a two-dimensional (2D) numpy array that contains x and y coordinates of
points in 2D space. Accordingly, this numpy array will have a dimension of number
of points × 2.
- **low th**: a scalar float quantity defines the lower accepted threshold values for the
Euclidean distance squared value. Use 0.95
- **up th**: a scalar float quantity defines the upper accepted threshold values for the
Euclidean distance squared value. Use 1.05

**Returns:**
Two numpy arrays, which are indexes of points where the distance squared is within the defined limit. the first numpy array contains the indexes for the first points and the second array contains the indexes for the secondset of points.

3. **plot_carbon_carbon**: The third function, plot carbon carbon, will plot two small dots (carbons) at a specific location and connected by a line. 

**Parameters:**
- p1 : numpy array coordinates of point 1
- p2 : numpy array coordinates of point 2

4. **plot_mesh**: Computes the actual mesh of the nanotube. This function will use plot carbon carbon to generate the figure with small dots, connecting every dot with a straight line.

**Parameters:**
- xx : numpy array indexes for the first points
- yy : numpy array indexes for the second points
- points : numpy array x and y coordinates of carbons in 2D space

## Instructions
1. Set up the graphene sheet orientation as described.
2. Use the functions in `nano_plot.py` to visually depict the graphene sheet and carbon nanotube.
3. Observe the relationships between atoms and their connectivity in the constructed nanotube.

## Conclusion
This project provided a hands-on approach to understanding the structure and formation of carbon nanotubes from graphene sheets, emphasizing the importance of orientation and connectivity in their formation.
