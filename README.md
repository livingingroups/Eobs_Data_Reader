
# Eobs_Data_Reader

The `Eobs_Data_Reader` function is an R-based tool designed to automate the pre-processing of **Acceleration**, **Magnetometry**, and **Orientation** data collected from Eobs devices. Leveraging `data.table` coding, this function accelerates the extraction and processing of large datasets, specifically designed to handle individual animal datasets downloaded directly from [Movebank](https://www.movebank.org/).

![Movebank Screenshot](Movebankimage.png)


**[This function is currently in development. Your feedback and testing are highly appreciated! If you encounter any issues or have suggestions for improvement, please feel free to share them.]**


## Key Features

### General Functionality
- **Subset Data**: Enables users to subset the dataset to a specified timestamp interval.
- **Data Conversion**: Converts row-wise, space-separated sensor readings (Acceleration, Magnetometry, Quaternions) into long format, where each row represents a single axis/component reading.
- **Dynamic Handling of Sampling Rates**: 
  - Provides an option to standardize acceleration data to a user-specified sampling frequency (Hz), ensuring consistency when multiple sampling rates exist in the dataset.
  - Handles cases where fewer than three acceleration axes are recorded.
- **Timestamp Extrapolation**: Automatically extrapolates timestamps within data bursts, ensuring uniformly spaced time intervals.

### Sensor-Specific Features
#### Acceleration Data:
- Automatically recognizes the sensor type (Legacy or IMU Accelerometer) and applies the appropriate transformation to convert raw acceleration data into standardized units of **g** (1 g ≈ 9.81 m/s²).
- Provides columns and visualizations for:
  - Acceleration burst durations and intervals over time.
  - Sensor type classification (Legacy vs IMU).
  - Burst IDs merged during continuous acceleration recordings.

- Computes per-burst acceleration metrics:
  - **Static & dynamic Acceleration**: Calculated using user-defined rolling averages.
  - **VeDBA**: Computes Vectorial Dynamic Body Acceleration for each burst, accounting for continuous sampling.

#### Orientation Data (Quaternions):
- Converts quaternions (if 20 Hz Orientation data is present) into mathematically meaningful components (W, X, Y, Z) for downstream analysis.

#### Magnetometry Data:
- Processes magnetometry data (if present) into long format.

### Advanced Error Handling and Progress Indicators
- Validates required columns for each data type, skipping sections with missing data while providing clear console messages.
- Flags potential duplicate timestamps in a dedicated `duplicate_times` column.
- Includes progress bars and detailed console messages to keep users informed about the current processing step.

### Flexible Outputs
  - Dynamically generates outputs based on the available data:
  - **Acceleration Data Only:** Acceleration Data Only: Returns a data frame named "Acceleration Data".
  - **All Data Types Available:** Acceleration Data Only: Returns a data frame named "Acceleration Data".
         - If Magnetometry Data and Quaternion Data timestamps match, these are combined into a data frame named "Orientation Data" and returned alongside the "Acceleration Data".
         - If timestamps do not match, a list is returned with relevant names: "Magnetometry Data", "Quaternion Data", and "Acceleration Data".
  - **Magnetometry and Quaternion Data Only:**
          - If timestamps match, the two datasets are combined into a single data frame named "Orientation Data".
          - If timestamps do not match, a list with "Magnetometry Data" and "Quaternion Data" is returned.
  - **Combination of Acceleration Data with Magnetometry or Quaternion Data:** Returns a list of the available datasets with appropriate names.

![Acceleration Summary Plots](Summaryplot.png)

## Example Workflow

1. **Download Data from Movebank**:
   - Ensure the dataset includes relevant columns for Acceleration, Magnetometry, and/or Orientation data.
   
2. **Load Data in R**:
   - This must be read in using **read.csv()** and **NOT** `fread()`
   ```r
   df <- read.csv("example_movebank_data.csv") 
   ```

3. **Run the Function**:
   ```r
   processed_data <- Eobs_Data_Reader(
     data = df,
     rolling_mean_width = 40, # If 20 Hz frequency, this corresponds to 2 s running mean
     standardised_freq_rate = 20,
     start_timestamp = "2024-05-01 00:00:00",
     end_timestamp = "2024-06-01 23:59:59",
     plot = TRUE
   )
   ```

### Required R Packages

The function relies on the following R packages. If any are not installed, the function will automatically install them before proceeding.

- **data.table** 
- **pbapply** 
- **tidyr** 
- **dplyr**
- **ggplot2**
- **cowplot**
- **viridis**

## License

This project is licensed under the MIT License.

## Contact

For questions, bug reports, suggestions, or contributions, please contact:
- Richard Gunner
- Email: rgunner@ab.mpg.de
- GitHub: [Richard6195](https://github.com/Richard6195)