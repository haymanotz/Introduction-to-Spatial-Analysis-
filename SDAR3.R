# Install packages if you don't have them already
 install.packages(c("sp","sf", "stars", "gstat", "spdep", "tmap"))

# Load necessary libraries
library(sp)      # For spatial data classes
library(gstat)   # The core geostatistics package (variograms, kriging)
library(spdep)   # For spatial autocorrelation (Moran's I, Geary's C)
library(tmap)    # For mapping results (optional, but nice!)

# Load the 'meuse' dataset
data(meuse)
# Convert the data.frame to a SpatialPointsDataFrame (required by 'sp' and 'gstat')
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("EPSG:28992") # Assign a coordinate system

OR
# Modern way to assign CRS in the 'sp' package
proj4string(meuse) <- CRS("EPSG:28992")

library(sf)
# If starting from a dataframe:
meuse_sf <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992)
meuse_sf

## 1. Create a Spatial Weights Matrix (e.g., k-nearest neighbors)
meuse.knn <- knearneigh(meuse, k=4)
meuse.nb <- knn2nb(meuse.knn)
meuse.listw <- nb2listw(meuse.nb)

# We'll analyze the 'zinc' concentration
zinc <- meuse$zinc

## 2. Moran's I (Test for **positive** spatial autocorrelation)
# H0: No spatial autocorrelation (random distribution)
moran_test <- moran.test(zinc, meuse.listw)
print(moran_test)

## 3. Geary's C (Test for **negative** spatial autocorrelation)
# H0: No spatial autocorrelation (random distribution)
geary_test <- geary.test(zinc, meuse.listw)
print(geary_test)

## 1. Calculate the Experimental Semi-Variogram
# 'log(zinc)' is often used because zinc concentration is highly skewed
vario.exp <- variogram(log(zinc) ~ 1, meuse) # ~ 1 indicates an ordinary variogram
plot(vario.exp)

## 2. Fit a Theoretical Variogram Model (e.g., Spherical)
# The 'vgm' function specifies the starting parameters
vario.model <- vgm(psill=0.5, model="Sph", range=1000, nugget=0.1)

# Fit the model to the experimental variogram
vario.fit <- fit.variogram(vario.exp, vario.model)
print(vario.fit)
plot(vario.exp, vario.fit) # Plot the data points and the fitted curve

# Calculate directional (anisotropic) variograms
vario.dir <- variogram(log(zinc) ~ 1, meuse, alpha=c(0, 45, 90, 135))
plot(vario.dir)

#4.1. Concept of Spatial Interpolation
#We need a grid to predict values onto. This grid is called the prediction map.
# Load the 'meuse.grid' dataset for the prediction locations
data(meuse.grid)
coordinates(meuse.grid) <- ~x+y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
gridded(meuse.grid) = TRUE # Convert to a SpatialPixelsDataFrame

# Interpolate 'log(zinc)' using IDW (power=2)
# formula: log(zinc) ~ 1 (simple IDW)
idw_zinc <- idw(log(zinc) ~ 1, meuse, meuse.grid, idp=2.0)

# Convert to a standard data frame for tmap plotting
idw_df <- as.data.frame(idw_zinc)
names(idw_df)[1] <- "log_zinc_IDW"

# Map the result
tm_shape(idw_zinc) +
  tm_raster(col="var1.pred", title="IDW log(Zinc)", style="quantile") +
  tm_dots(col="log(zinc)", size=0.1, data=meuse)


#OR
library(stars)

# 1. Convert the IDW result to a stars object
idw_stars <- st_as_stars(idw_zinc)

# Ensure the column exists in the spatial object
meuse_sf$log_zinc <- log(meuse_sf$zinc)

# Plotting with tmap v4
tm_shape(idw_stars) + 
  tm_raster(
    col = "var1.pred", 
    col.scale = tm_scale_intervals(style = "quantile", values = "viridis"),
    col.legend = tm_legend(title = "IDW log(Zinc)")
  ) +
tm_shape(meuse_sf) + 
  tm_dots(
    fill = "log_zinc",  # Now this column is guaranteed to exist
    size = 0.1,
    fill.scale = tm_scale_continuous(values = "magma"), # Changed palette so dots are visible
    fill.legend = tm_legend(title = "Actual Samples")
  )
# Use the fitted variogram model from 3.4
# vario.fit

# Interpolate using Ordinary Kriging
# formula: log(zinc) ~ 1 (Ordinary Kriging uses a constant mean)
ok_zinc <- krige(log(zinc) ~ 1, meuse, meuse.grid, model=vario.fit)

# Convert to a standard data frame for tmap plotting
ok_df <- as.data.frame(ok_zinc)
names(ok_df)[1] <- "log_zinc_OK_Prediction"

# Map the Kriging prediction result
tm_shape(ok_zinc) +
  tm_raster(col="var1.pred", title="OK log(Zinc) Prediction", style="quantile") +
  tm_dots(col="log(zinc)", size=0.1, data=meuse)

#OR 
library(stars)
library(sf)

# 1. Convert Kriging results to stars (for the raster)
ok_stars <- st_as_stars(ok_zinc)

# 2. Ensure your points are in sf format and have the log variable
meuse_sf <- st_as_sf(meuse)
meuse_sf$log_zinc <- log(meuse_sf$zinc)

# 3. Plot the Prediction Map
tm_shape(ok_stars) + 
  tm_raster(
    col = "var1.pred", 
    col.scale = tm_scale_intervals(style = "quantile", values = "viridis"),
    col.legend = tm_legend(title = "OK log(Zinc) Prediction")
  ) +
tm_shape(meuse_sf) + 
  tm_dots(
    fill = "log_zinc", 
    size = 0.1, 
    fill.scale = tm_scale_continuous(values = "magma"),
    fill.legend = tm_legend(title = "Actual Samples")
  )

#4.3. Uncertainty Assessment in Spatial Interpolation
#A major advantage of Kriging is that it simultaneously provides the Kriging Variance (uncertainty/error map).
# The Kriging function ('krige') automatically calculates the prediction variance
ok_variance <- ok_zinc
names(ok_variance)[2] <- "log_zinc_OK_Variance"

# Map the Kriging Variance (uncertainty)
tm_shape(ok_variance) +
  tm_raster(col="var1.pred", title="OK Prediction Variance", style="cont")
  # High variance (darker areas) indicates higher uncertainty, often far from sample points


# 1. Check the actual column names in your variance object
names(ok_variance)


library(stars)

# 1. Convert the result to a stars object if it isn't one
ok_variance_stars <- st_as_stars(ok_variance)

# 2. Map the Kriging Variance
tm_shape(ok_variance_stars) +
  tm_raster(
    col = "var1.pred", # Changed from .pred to .var for variance/uncertainty
    col.scale = tm_scale_continuous(values = "brewer.yl_or_rd"), # Updated palette name for v4
    col.legend = tm_legend(title = "OK Prediction Variance")
  )

# Note: 'meuse.grid' also contains the 'dist' variable
uk_zinc <- krige(log(zinc) ~ dist, meuse, meuse.grid, model=vario.fit)

# Map the UK prediction
tm_shape(uk_zinc) +
  tm_raster(col="var1.pred", title="UK log(Zinc) Prediction (~dist)", style="quantile")

Or
library(stars)

# 1. Convert the Universal Kriging result to stars
uk_stars <- st_as_stars(uk_zinc)

# 2. Map using tmap v4 syntax
tm_shape(uk_stars) +
  tm_raster(
    col = "var1.pred", 
    col.scale = tm_scale_intervals(style = "quantile", values = "YlOrBr"),
    col.legend = tm_legend(title = "UK log(Zinc) (~dist)")
  ) +
  tm_layout(main.title = "Universal Kriging Predictions", frame = FALSE)




#Geostatistical Analysis Code (Using Your Data)
#Setup and Data Loading
#This section ensures your CSV data is loaded correctly and converted into the necessary spatial format (SpatialPointsDataFrame).
# 1. Install and Load Packages (Run only if needed)
install.packages(c("sp", "gstat", "spdep", "tmap"))

library(sp)      # For spatial data classes
library(gstat)   # Core geostatistics package
library(spdep)   # For spatial autocorrelation tests
library(tmap)    # For mapping results (optional)
Read Your Data
# Read CSV file
data <- read.csv("C:/Users/Admin/Desktop/sppd.csv")

# View structure
str(data)
head(data)
#Convert CSV to Spatial Data
coordinates(data) <- ~longitude+latitude
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")

plot(data, pch=16, col="blue", main="Spatial Distribution of Data")

#3.2 Spatial Autocorrelation
#Moran’s I
# Create spatial neighbors (distance-based)
nb <- dnearneigh(coordinates(data), 0, 2)   # adjust distance if needed
lw <- nb2listw(nb, style="W")

# Moran's I
moran.test(data$value, lw)

Geary’s C
geary.test(data$value, lw)

3.3 Isotropy vs Anisotropy
# Directional variograms
vgm_dir <- variogram(popDensity ~ 1, data, alpha=c(0,45,90,135))

plot(vgm_dir, main="Directional Variograms (Anisotropy Test)")

csv_file_path <- "C:/Users/Admin/Desktop/earthquake_data_tsunami.csv" 
X_COORD_COL <- "longitude" 
Y_COORD_COL <- "latitude" 
CRS_CODE    <- 4326        # EPSG code for WGS84 (Latitude/Longitude)
Z_VALUE_COL <- "magnitude"       # The column for the raster value
library(readr) 
library(sf)
library(terra)

# 1. Define File Path and Coordinates (ADJUST THESE)
csv_file_path <- "C:/Users/Admin/Desktop/earthquake_data_tsunami.csv" 
X_COORD_COL <- "longitude" 
Y_COORD_COL <- "latitude" 
CRS_CODE    <- 4326        # EPSG code for WGS84 (Latitude/Longitude)
Z_VALUE_COL <- "magnitude"       # The column for the raster value
# 2. Read the CSV and convert it to a spatial point object (sf)
point_data_df <- read_csv(csv_file_path)
# This creates the VECTOR data model (sf object)
vector_points_sf <- st_as_sf(
  point_data_df,
  coords = c(X_COORD_COL, Y_COORD_COL),
  crs = CRS_CODE
)
proj4string(vector_points_sf) <- CRS("+proj=longlat +datum=WGS84")

plot(vector_points_sf, pch=16, col="blue", main="Spatial Distribution of Data")
# This plots ONLY the points, avoiding the 10-attribute warning
plot(st_geometry(vector_points_sf), pch=16, col="blue", main="Spatial Distribution of Data")


# 1. Define Spatial Weights (e.g., 4-Nearest Neighbors)
data.knn <- knearneigh(data, k=4)
data.nb <- knn2nb(data.knn)
data.listw <- nb2listw(data.nb)

# 2. Extract the Variable of Interest
# Note: Use get(Z_Variable_Col) to dynamically reference the column name
Z_variable <- data[[Z_Variable_Col]] 

## Moran's I Test
moran_test <- moran.test(Z_variable, data.listw)
print(moran_test)

## Geary's C Test
geary_test <- geary.test(Z_variable, data.listw)
print(geary_test)

# Directional variograms
vgm_dir <- variogram(Z_variable ~ 1, data, alpha=c(0,45,90,135))

plot(vgm_dir, main="Directional Variograms (Anisotropy Test)")

#Empirical Variogram
vgm_emp <- variogram(Z_variable ~ 1, data)
plot(vgm_emp, main="Empirical Variogram")

Fit Variogram Model
vgm_model <- fit.variogram(vgm_emp, model=vgm("Sph"))

plot(vgm_emp, vgm_model, main="Fitted Variogram Model")


#Create Prediction Grid
# Create grid
grd <- expand.grid(
  x = seq(min(data@coords[,1]), max(data@coords[,1]), length=100),
  y = seq(min(data@coords[,2]), max(data@coords[,2]), length=100)
)

coordinates(grd) <- ~x+y
gridded(grd) <- TRUE
proj4string(grd) <- CRS("+proj=longlat +datum=WGS84")

#Inverse Distance Weighting (IDW)
idw_pred <- idw(Z_variable ~ 1, data, grd)

spplot(idw_pred["var1.pred"], main="IDW Interpolation")
# Check for duplicates
any(duplicated(coordinates(data)))

# If TRUE, remove or aggregate:
data <- zerodist(data, zero = 0.0) # removes duplicates
# 1. Convert the matrix back to a data frame
data_df <- as.data.frame(data)

# 2. Tell R which columns are the coordinates (replace with your actual column names)
# If your columns are named 'x' and 'y', use:
coordinates(data_df) <- ~V1+V2
 
colnames(data_df)
# 3. Now try the variogram again
vario_emp <- variogram(Z_variable ~ 1, data = data_df)
plot(vario_emp, model = vgm_model)
# Now plot the empirical points and your model together
plot(vario_emp, model = vgm_model)

# Plot to see if the blue line actually follows the points
plot(variogram(Z_variable ~ 1, data), vgm_model)

#Kriging Interpolation

kriging_pred <- krige(Z_variable ~ 1, data, grd, model=vgm_model)

spplot(kriging_pred["var1.pred"], main="Kriging Prediction")

4.3 Uncertainty Assessment
spplot(kriging_pred["var1.var"], main="Kriging Prediction Variance")

High variance = high uncertainty