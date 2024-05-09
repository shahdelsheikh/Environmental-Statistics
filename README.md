# Environmental Statistics (Spatial Data Analysis for Temperature in California) using R.
## Introduction
One of the most substantial environmental factorsin the natural world is temperature, affecting the whole ecosystem with all living things. Understanding the patterns and variations of environmental phenomena such as temperature across geographical regions relies heavily on spatial analysis using 456 locations in California

## Objectives
The main objectives of this project are as follows:

1. Analyze the average temperature patterns in California during the second quarter of the year.
2. Explore the spatial relationships, distributions, and trends of temperature across geographical regions.
3. Identify spatial autocorrelation patterns to understand the dependency structure and clustering of temperature values.
4. Generate interpolated maps to assess precision levels of temperature estimation in different areas of California.

## Methodology
### Data exploration 
It included graphical methods (Box and Bubble plot) and the geographical weighted regression model (GWRM) as a method of exploration. GWR is a non- stationary technique that models spatially varying relationships. 
### Spatial Autocorrelation 
Global Moran's I and Localized Moran to assess how similar or dissimilar spatial data are in terms of average temperature across different regions of California based on their geographic proximity.
### Spatial Interpolation
#### Inverse Distance Weighting
IDW was utilizied to estimate the average temperature values at unsampled locations. The optimal
value for power parameter was equal 3.881, knowing that we rarely choose p greater than 5, the optimal value of p is large in comparison to 5, indicating that the spatial dependency structure in the second quarter is local.
#### Kriging 
Kriging was utilized to estimate the average temperature values across California. The Exponential model was selected based on its ability to capture spatial correlation patterns after constructing variograms and finding out that it had the least MSE value. Then, we created interpolated maps to determine where to locate our next stations.
