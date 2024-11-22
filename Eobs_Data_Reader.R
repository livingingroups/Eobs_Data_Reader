# Eobs_Data_Reader

## Function Description ##
# The Eobs_Data_Reader function is an R script designed to preprocess high-frequency Acceleration, 
# Magnetometry, and Orientation data collected from Eobs devices. It streamlines data handling, and cleaning, 
# for datasets downloaded directly from Movebank, allowing users to focus on analysis rather 
# than data preparation.

## Function Inputs ##
 # df (Required): A data.frame containing raw sensor data downloaded from Movebank, with original column names.
 # rolling_mean_width (Default: 40): The window size (in samples) used to compute rolling averages for static acceleration.
 # standardised_freq_rate (Default: NULL): A target frequency (in Hz) to which acceleration data should be downsampled.
 # start_timestamp (Default: NULL): The earliest timestamp to include in the dataset (in YYYY-MM-DD HH:MM:SS.sss format).
 # end_timestamp (Default: NULL): The latest timestamp to include in the dataset (in YYYY-MM-DD HH:MM:SS.sss format).
 # plot (Default: TRUE): A logical flag indicating whether to generate visualizations of the processed data. If TRUE, the function produces plots summarizing acceleration burst durations, sampling intervals, and sensor type classifications.
 
# Note: While fread() is faster than read.csv(), it may alter column names in unexpected ways.
# To ensure the function correctly interprets the column names, use read.csv() instead.
# Example:  df <- read.csv("C:/Users/richa/OneDrive/Desktop/Eobs data Converted/Albus/Albus (1).csv")



Eobs_Data_Reader <- function(df, rolling_mean_width = 40, standardised_freq_rate = NULL, start_timestamp = NULL, end_timestamp = NULL, plot = TRUE){
  
  suppressWarnings ({
    # Suppress dplyr summarise info
    options(dplyr.summarise.inform = FALSE)
    
    ### Automatically Load or Install Required Packages ###
    required_packages <- c("dplyr", "tidyr", "ggplot2", "data.table", "pbapply", "cowplot", "viridis")
    
    install_and_load <- function(packages) {
      for (pkg in packages) {
        if (!require(pkg, character.only = TRUE)) {
          message(paste("Installing missing package:", pkg))
          install.packages(pkg, dependencies = TRUE, type = "source")
          suppressMessages(library(pkg, character.only = TRUE))
        } else {
          suppressMessages(library(pkg, character.only = TRUE))
        }
      }
    }
    missing_pkgs <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
    if (length(missing_pkgs) > 0) {
      message(paste("Installing missing packages:", paste(missing_pkgs, collapse = ", ")))
    }
    # Execute the function to load all required packages
    install_and_load(required_packages)
    
    # Function to expand raw IMU space-seperated data into long format
    process_data <- function(params) {
      row_id <- params$row_id
      raw_values <- params$raw_values
      axes <- params$axes
      if (is.na(raw_values) || !is.character(raw_values)) {
        return(data.table(row_id = row_id, X = NA, Y = NA, Z = NA))
      }
      # Split raw values into a numeric vector
      values <- as.numeric(unlist(strsplit(raw_values, " ")))
      # Determine axis order based on `axes`
      axis_order <- switch(
        axes,
        "XYZ" = c("X", "Y", "Z"),
        "XY" = c("X", "Y"),
        "XZ" = c("X", "Z"),
        "YZ" = c("Y", "Z"),
        "X" = c("X"),
        "Y" = c("Y"),
        "Z" = c("Z"),
        stop("Unknown axes configuration")
      )
      # Repeat axis order to match the length of the raw values
      repeat_order <- rep(axis_order, length.out = length(values))
      # Create a data.table with row_id and expanded acceleration data
      data.table(
        row_id = row_id,
        X = if ("X" %in% repeat_order) values[repeat_order == "X"] else NA,
        Y = if ("Y" %in% repeat_order) values[repeat_order == "Y"] else NA,
        Z = if ("Z" %in% repeat_order) values[repeat_order == "Z"] else NA
      )
    }
    
    # Custom function
    downsample_data_table <- function(data, target_rate) {
      # Validate input
      if (!"eobs.acceleration.sampling.frequency.per.axis" %in% names(data)) {
        stop("Column 'eobs.acceleration.sampling.frequency.per.axis' is missing in the data.")
      }
      # Downsample each burst independently
      downsampled_data <- data[, {
        original_rate <- unique(eobs.acceleration.sampling.frequency.per.axis)
        if (length(original_rate) != 1) {
          stop("Multiple sampling frequencies found within the same burst. Please check data integrity.")
        }
        if (target_rate > original_rate) {
          stop("standardised_freq_rate must be less than or equal to the lowest original sampling frequency.")
        }
        # Calculate downsampling factor
        factor <- original_rate / target_rate
        # Select indices based on downsampling factor
        indices <- seq(1, .N, by = factor)
        rounded_indices <- round(indices)
        valid_indices <- rounded_indices[rounded_indices > 0 & rounded_indices <= .N]
        # Return only the downsampled rows
        .SD[valid_indices]
      }, by = row_id]
      
      return(downsampled_data)
    }
    
    # Helper function to check required columns
    validate_columns <- function(data, required_columns, section_name) {
      if (is.null(required_columns) || length(required_columns) == 0) {
        stop(paste("No columns specified for validation in", section_name, "section."))
      }
      missing_columns <- setdiff(required_columns, colnames(data))
      if (length(missing_columns) > 0) {
        message(paste("Skipping", section_name, "section. Missing columns:", paste(missing_columns, collapse = ", ")))
        return(FALSE)
      }
      return(TRUE)
    }
    
    # Validate initial required columns
    initial_required_columns <- c(
      "timestamp", 
      "individual.taxon.canonical.name", 
      "tag.local.identifier", 
      "individual.local.identifier",
      "sensor.type"
    )
    
    ########################################################################################################################
    # Convert to data.table for efficient processing
    setDT(df) 
    
    if (!validate_columns(df, initial_required_columns, "initial formatting/subsetting")) {
      stop("Data is missing essential columns for initial formatting. Cannot proceed.")
    }
    
    # Ensure timestamps are properly formatted
    message("Formatting/Subsetting timestamps...")
    options(digits.secs = 3)  # Ensure decimal seconds are displayed
    # Convert 'timestamp' column
    df[, timestamp := as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS")]
    
    # Validate start and end timestamps
    if (!is.null(start_timestamp)) {
      start_timestamp <- as.POSIXct(start_timestamp, format = "%Y-%m-%d %H:%M:%OS")
      if (is.na(start_timestamp)) stop("Invalid start_timestamp format. Please use 'YYYY-MM-DD HH:MM:SS.sss'.")
    }
    
    if (!is.null(end_timestamp)) {
      end_timestamp <- as.POSIXct(end_timestamp, format = "%Y-%m-%d %H:%M:%OS")
      if (is.na(end_timestamp)) stop("Invalid end_timestamp format. Please use 'YYYY-MM-DD HH:MM:SS.sss'.")
    }
    
    # Subset data between start and end timestamps, if provided
    if (!is.null(start_timestamp) || !is.null(end_timestamp)) {
      message("Subsetting data based on provided start and/or end timestamps...")
      df <- df[
        (!is.null(start_timestamp) & timestamp >= start_timestamp) &
          (!is.null(end_timestamp) & timestamp <= end_timestamp)
      ]
      if (nrow(df) == 0) stop("No data remains after subsetting with the given timestamps.")
    }
    
    ################################################################################################################################################################################################
    ################################################################################################################################################################################################
    #################################################ACCELERATION###################################################################################################################################
    
    # Check if "acceleration" exists in the sensor.type column
    if ("acceleration" %in% df$sensor.type) {
      required_columns <- c(
        "timestamp", 
        "individual.taxon.canonical.name", 
        "tag.local.identifier", 
        "individual.local.identifier",
        "eobs.acceleration.axes", 
        "eobs.acceleration.sampling.frequency.per.axis", 
        "eobs.accelerations.raw"
      )
      if (validate_columns(df, required_columns, "Acceleration")) {
        # Check if there are any non-NA space-separated character strings in eobs.accelerations.raw
        if (any(!is.na(df$eobs.accelerations.raw) & nzchar(df$eobs.accelerations.raw) & !is.na(df$eobs.acceleration.sampling.frequency.per.axis) & !is.na(df$eobs.acceleration.axes))) {
          # Define columns required
          
          message("Starting the acceleration data processing pipeline...")
          # Filter rows for acceleration data
          acc <- df[sensor.type == "acceleration"]
          
          # Arrange rows by timestamp
          acc <- acc[order(timestamp)]
          rownames(acc) <- NULL
          
          # Select and Reorder
          acc <- acc[, required_columns[required_columns %in% names(acc)], with = FALSE]
          setcolorder(acc, required_columns[required_columns %in% names(acc)])
          
          message(paste("Extracting acceleration to long format..."))
          
          # Prepare parameters for parallel processing
          params_list <- lapply(seq_len(nrow(acc)), function(i) {
            list(
              row_id = i,
              raw_values = acc$eobs.accelerations.raw[i],
              axes = acc$eobs.acceleration.axes[i]
            )
          })
          
          # Process acceleration data to long format
          processed_acc <- pblapply(params_list, process_data)
          # Combine results into a single data.table
          processed_acc <- rbindlist(processed_acc, fill = TRUE)
          
          # Add metadata back using the row_id
          processed_acc <- merge(processed_acc, acc[, .(row_id = .I, timestamp,
                                                        individual.taxon.canonical.name, tag.local.identifier, 
                                                        individual.local.identifier, eobs.acceleration.axes, 
                                                        eobs.acceleration.sampling.frequency.per.axis)], 
                                 by = "row_id", all.x = TRUE)
          
          # Remove redundant acc data.table
          rm(acc)
          # Identify and dynamically handle present axes
          present_axes <- intersect(c("X", "Y", "Z"), names(processed_acc))
          setcolorder(processed_acc, c(setdiff(names(processed_acc), present_axes), present_axes)) # All other columns except the present axes
          # Rename the present axes to `acc_x`, `acc_y`, and `acc_z`
          new_axis_names <- paste0("acc_", tolower(present_axes))  # Convert to acc_x, acc_y, acc_z
          setnames(processed_acc, old = present_axes, new = new_axis_names)
          
          # Assign sensor type based on burst length
          message("Assigning burst lengths...")
          processed_acc[, burst_length := .N, by = row_id]
          # Add 'sensor_type' column based on burst length
          processed_acc[, sensor_type := ifelse(burst_length == 24, "IMU", "Legacy")] # One IMU ACC dataset consists of 24*3 samples
          
          # Downsample the processed acceleration data ?
          if (!is.null(standardised_freq_rate)) {
            message("Downsampling data to target rate: ", standardised_freq_rate, " Hz")
            processed_acc <- downsample_data_table(processed_acc, target_rate = standardised_freq_rate)
            processed_acc[, eobs.acceleration.sampling.frequency.per.axis := standardised_freq_rate] # New sampling rate
            processed_acc[, burst_length := .N, by = row_id] # Redo burst length
          }
          
          # Identify which axes are present in the processed data
          present_axes <- intersect(c("acc_x", "acc_y", "acc_z"), names(processed_acc))
          
          ###############################
          message("Extrapolating within burst timestamps...")
          # Extrapolate timestamps within each burst using row_id
          processed_acc[, extrapolated_timestamp := timestamp]  # Start with the original timestamp
          
          # Add time increment based on the sampling frequency
          processed_acc[, time_increment := 1 / eobs.acceleration.sampling.frequency.per.axis, by = row_id]
          
          # Get the first timestamp and sensor type of each burst and the next burst's details
          burst_starts <- processed_acc[, .(burst_start = first(timestamp), sensor_type = first(sensor_type)), by = row_id]
          burst_starts[, next_burst_start := shift(burst_start, type = "lead")]
          burst_starts[, next_sensor_type := shift(sensor_type, type = "lead")]
          # Remove sensor_type
          burst_starts[, c("sensor_type") := NULL]
          
          # Merge the burst start times and sensor types back into the main dataset
          processed_acc <- merge(processed_acc, burst_starts, by = "row_id", all.x = TRUE)
          
          # Initialize progress counter
          total_rows <- length(unique(processed_acc$row_id))
          progress_step <- floor(total_rows / 100)  # Update every 1% of progress
          message("Adjusting extrapolated timestamps...")
          
          processed_acc[, extrapolated_timestamp := {
            # Use a static counter to track progress
            if (.GRP %% progress_step == 0) {
              message(paste0("Progress: ", round(.GRP / total_rows * 100), "%"))
            }
            
            if (length(unique(sensor_type)) == 1) {
              # Handle single sensor type dataset
              timestamp + (seq_len(.N) - 1) * time_increment
            } else {
              {
                # Check if the next burst starts within 1.2 seconds and the sensor type remains the same
                if (!is.na(next_burst_start[1]) &&
                    (next_burst_start[1] - burst_start[1] <= 1.2) &&
                    sensor_type[1] == next_sensor_type[1]) {
                  # Linearly interpolate timestamps up to just before the next burst's start
                  seq(
                    from = burst_start[1], 
                    to = next_burst_start[1] - (1 / eobs.acceleration.sampling.frequency.per.axis[1]), 
                    length.out = .N
                  )
                } else {
                  # Default extrapolation behavior
                  timestamp + (seq_len(.N) - 1) * time_increment
                }
              }
            }
          }, by = row_id]
          
          # Drop intermediate columns if necessary
          processed_acc[, c("time_increment", "burst_start", "next_burst_start", "next_sensor_type") := NULL]
          
          # Add a time difference column (not grouped by row_id)
          processed_acc[, time_diff := c(NA, diff(extrapolated_timestamp))]
          
          # Mark duplicates in the `extrapolated_timestamp` column
          processed_acc[, duplicate_times := duplicated(extrapolated_timestamp)]
          
          ################# Plot uninterrupted duration of each sensor type ########################
          message("Processing sampling durations and uniterupted burst lengths...")
          # Subset data to include only the first row of each burst
          burst_data <- processed_acc[, .SD[1], by = row_id]  # Keep the first row of each burst
          
          # Add a group identifier for uninterrupted sequences of the same sensor type
          burst_data[, sensor_sequence := rleid(sensor_type)]  # Create unique ID for runs of the same sensor type
          
          # Compute sampling intervals only within uninterrupted sequences
          burst_data[, sampling_interval := c(NA, diff(as.numeric(timestamp))), by = sensor_sequence]  # Compute intervals
          
          # Merge sampling_interval back into processed_acc
          processed_acc <- merge(
            processed_acc, 
            burst_data[, .(row_id, sensor_sequence, sampling_interval)],  # Select only relevant columns from burst_data
            by = "row_id", 
            all.x = TRUE  # Preserve all rows in processed_acc
          )
          
          # Initialize burst_id
          processed_acc[, burst_id := {
            indices <- seq_len(.N)
            # Initialize burst_id
            burst_id <- 1
            # Create a vector to hold burst IDs for the group
            burst_ids <- integer(.N)
            burst_ids[1] <- burst_id  # Assign the first burst ID
            
            # Iterate through rows to assign burst_id based on conditions
            pbapply::pblapply(indices[-1], function(i) {
              # Assign burst_id based on conditions
              if (is.na(sensor_sequence[i]) || is.na(sensor_sequence[i - 1])) {
                burst_ids[i] <<- burst_id  # Retain the same burst_id
              } else if (time_diff[i] > 0.5 || sensor_sequence[i] != sensor_sequence[i - 1]) { #0.5 s the cut off for continuous bursts
                burst_id <<- burst_id + 1  # Increment burst_id when conditions are met
                burst_ids[i] <<- burst_id
              } else {
                burst_ids[i] <<- burst_id  # Keep the same burst_id
              }
            })
            
            burst_ids  # Return the vector of burst IDs
          }]
          
          if (plot == TRUE) {
            ########################### Summary plots #################################
            
            message("Summary plots of sensor type classification and sampling intervals")
            
            # Dynamically compute features for each burst based on present axes
            burst_features <- processed_acc[, c(
              # Compute mean acceleration only for present axes
              lapply(present_axes, function(axis) list(mean = mean(get(axis), na.rm = TRUE))),
              # Add burst length and sampling frequency
              list(
                burst_length = .N,
                sampling_frequency = unique(eobs.acceleration.sampling.frequency.per.axis),
                sensor_type = unique(sensor_type)
              )
            ), by = row_id]
            
            # Identify indices of unnamed or default columns
            blank_col_indices <- which(grepl("^\\.|^$", colnames(burst_features)))
            
            # Ensure that the number of detected unnamed columns matches `present_axes`
            if (length(blank_col_indices) == length(present_axes)) {
              # Rename the unnamed columns dynamically
              colnames(burst_features)[blank_col_indices] <- paste0("mean_", present_axes)
            } else {
              # Throw a detailed error if there's a mismatch
              stop("Mismatch: Number of unnamed columns does not match the number of detected axes.\n",
                   "Unnamed column indices: ", paste(blank_col_indices, collapse = ", "), "\n",
                   "Detected axes: ", paste(present_axes, collapse = ", "))
            }
            
            # Ensure all measure.vars columns are numeric
            for (col in paste0("mean_", present_axes)) {
              burst_features[, (col) := as.numeric(get(col))]
            }
            
            # Melt data for faceted plotting
            plot_data <- melt(
              burst_features,
              id.vars = c("burst_length", "sampling_frequency", "sensor_type"),
              measure.vars = paste0("mean_", present_axes),
              variable.name = "axis",
              value.name = "mean_acc"
            )
            
            #Convert back to data.table
            setDT(plot_data)
            
            # Compute burst durations
            burst_summary <- processed_acc[, .(
              burst_start_time = min(extrapolated_timestamp, na.rm = TRUE),
              burst_end_time = max(extrapolated_timestamp, na.rm = TRUE),
              burst_duration = max(extrapolated_timestamp, na.rm = TRUE) - min(extrapolated_timestamp, na.rm = TRUE),
              sensor_type = first(sensor_type)
            ), by = burst_id]
            
            # Convert burst_duration to numeric (in seconds)
            burst_summary[, burst_duration := as.numeric(burst_duration)]
            
            # Define axis labels dynamically based on present axes
            axis_labels <- setNames(
              paste0("Mean Acceleration (", toupper(sub("mean_acc_", "", paste0("mean_", present_axes))), ")"),
              paste0("mean_", present_axes)
            )
            
            # Plot mean acceleration over sampling duration with dynamic axes
            p0 <- ggplot(plot_data, aes(x = burst_length / sampling_frequency, y = mean_acc, color = sensor_type)) +
              geom_boxplot() +
              facet_wrap(~axis, scales = "free_y", labeller = labeller(axis = axis_labels)) +  # Use dynamic axis labels
              labs(
                title = "Sensor Type Classification",
                x = "Burst Length (s)",
                y = "Mean (Analogue Digital) \nAcceleration",
                color = "Sensor Type"
              ) +
              scale_color_manual(values = c("Legacy" = "blue", "IMU" = "red")) +
              theme_minimal()
            
            # Plot sampling intervals
            p1 <- ggplot(burst_data, aes(x = extrapolated_timestamp, y = sampling_interval, colour = sensor_type)) +
              geom_point(alpha = 0.4) +
              labs(
                title = "Sampling Intervals Within Uninterrupted Sensor Sequences",
                x = "Timestamp",
                y = "Sampling Interval \n(s)",
                colour = "Sensor Type"
              ) +
              scale_color_manual(values = c("Legacy" = "blue", "IMU" = "red")) +
              theme_minimal()
            
            # Aggregate sampling intervals by half-hour and sensor type
            burst_data[, half_hour := as.POSIXct(
              cut(timestamp, breaks = "30 min"), , format = "%Y-%m-%d %H:%M:%OS"
            )]
            
            facet_data <- burst_data[, .(
              avg_sampling_interval = mean(sampling_interval, na.rm = TRUE)
            ), by = .(half_hour, sensor_type)]
            
            # Remove NA sampling intervals
            facet_data <- facet_data[!is.na(avg_sampling_interval)]
            
            # Create the heatmap-like plot
            p2 <- ggplot(facet_data, aes(x = half_hour, y = sensor_type, fill = avg_sampling_interval)) +
              geom_tile() +
              scale_fill_viridis_c(
                name = "Sampling Interval (s)",
                option = "plasma",
                na.value = "grey50"
              ) +
              labs(
                title = "Sampling Intervals by Sensor Type and Time",
                x = "Time (Half-Hour Bins)",
                y = "Sensor Type"
              ) +
              theme_minimal() +
              theme(
                plot.background = element_rect(fill = "black", color = NA),  # Black background
                panel.background = element_rect(fill = "black", color = NA),  # Black panel background
                axis.text = element_text(color = "white"),  # White axis text
                axis.title = element_text(color = "white"),  # White axis titles
                plot.title = element_text(color = "white", size = 14, face = "bold"),  # White plot title
                legend.background = element_rect(fill = "black", color = NA),  # Black legend background
                legend.text = element_text(color = "white", size = 8),  # White legend text
                legend.title = element_text(color = "white", size = 10),  # White legend title
                legend.key.height = unit(0.3, "cm"),  # Reduce height of legend scale bar
                legend.key.width = unit(0.3, "cm"),  # Reduce width of legend scale bar
                axis.text.x = element_text(angle = 45, hjust = 1, color = "white")  # Angled x-axis text in white
              )
            
            # Plot burst durations over time
            p3 <- ggplot(burst_summary, aes(x = burst_start_time, y = burst_duration, color = sensor_type)) +
              geom_point(size = 2, alpha = 0.7) +
              geom_line(aes(group = sensor_type), alpha = 0.5) +  # Optionally connect points by sensor type
              labs(
                title = "Burst Durations Over Time by Sensor Type",
                x = "Time",
                y = "Burst Duration (s)",
                color = "Sensor Type"
              ) +
              scale_color_manual(values = c("Legacy" = "blue", "IMU" = "red")) +
              theme_minimal() +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major = element_line(color = "grey90"),
                panel.grid.minor = element_blank()
              ) +
              facet_wrap(.~sensor_type, scales =  "free_y")
            
            
            # Combine the plots vertically
            combined_plot <- plot_grid(
              p0,                     # First plot (mean acc and burst duration)
              p1,                     # Second plot (sampling intervals)
              p2,                     # Third plot (heat map of sampling intervals)
              p3,                     # Fourth plot (continuous burst durations over time)
              ncol = 1,               # Arrange in one column (vertical stacking)
              align = "v",            # Align vertically
              rel_heights = c(1, 1, 1, 1) # Adjust relative heights if needed
            )
          }
          if (plot == TRUE) {
            # Explicitly print the combined plot
            print(combined_plot)
            rm(burst_summary, plot_data, burst_features, blank_col_indices, axis_labels, facet_data)
          }
          
          ##################################################################################################
          
          message("Transforming acceleration values into units of g and computing VeDBA...")
          
          # Transform acceleration to 'g' units for present axes only
          processed_acc[, (paste0(present_axes, "_g")) := lapply(present_axes, function(axis) {
            ifelse(sensor_type == "Legacy",
                   (get(axis) - 2048) / 512,
                   (get(axis) - 2048) / 1024)
          })]
          
          # Compute static and dynamic acceleration and VeDBA per burst_id using `frollmean`
          processed_acc[, c(
            paste0(present_axes, "_static"),
            paste0(present_axes, "_dynamic"),
            "VeDBA"
          ) := {
            # Compute static acceleration (rolling mean) for present axes using `frollmean`
            static <- lapply(present_axes, function(axis_g) {
              frollmean(get(paste0(axis_g, "_g")), n = rolling_mean_width, align = "center", na.rm = FALSE)
            })
            
            # Compute dynamic acceleration (subtract static from raw)
            dynamic <- Map(function(raw, static) raw - static, 
                           lapply(present_axes, function(axis_g) get(paste0(axis_g, "_g"))), 
                           static)
            
            # Compute VeDBA (Vectorial Dynamic Body Acceleration)
            vedba <- sqrt(Reduce(`+`, lapply(dynamic, function(d) ifelse(is.na(d), NA, d^2))))
            
            c(static, dynamic, list(vedba))
          }, by = burst_id]
          
          # Remove redundant variables
          rm(burst_features, present_axes, new_axis_names, burst_starts, burst_data)

          gc()
          
        } 
      }
    }
    ################################################################################################################################################################################################
    ################################################################################################################################################################################################
    #################################################MAGNETOMETERY##################################################################################################################################
    
    # Check if "magnetometer" exists in the sensor.type column
    if ("magnetometer" %in% df$sensor.type) {
      required_columns <- c(
        "timestamp", 
        "individual.taxon.canonical.name", 
        "tag.local.identifier", 
        "individual.local.identifier",
        "mag.magnetic.field.sampling.frequency.per.axis", 
        "mag.magnetic.fields.raw"
      )
      if (validate_columns(df, required_columns, "Magnetometer")) {
        # Check if there are any non-NA space-separated character strings in mag.magnetic.fields.raw
        if (any(!is.na(df$mag.magnetic.fields.raw) & nzchar(df$mag.magnetic.fields.raw) & !is.na(df$mag.magnetic.field.sampling.frequency.per.axis))) {
          
          message("Starting the magnetometery data processing pipeline...")
          
          # Filter rows for magnetometer data
          mag <- df[sensor.type == "magnetometer"]
          # Arrange rows by timestamp
          mag <- mag[order(timestamp)]
          rownames(mag) <- NULL
          # Make an XYZ column
          axes = "XYZ"
          mag[, eobs.magnetometery.axes := axes] 
          
          # Select and reorder
          required_columns <- c(
            "timestamp", 
            "individual.taxon.canonical.name", 
            "tag.local.identifier", 
            "individual.local.identifier",
            "eobs.magnetometery.axes",
            "mag.magnetic.field.sampling.frequency.per.axis", 
            "mag.magnetic.fields.raw"
          )
          mag <- mag[, required_columns[required_columns %in% names(mag)], with = FALSE]
          setcolorder(mag, required_columns[required_columns %in% names(mag)])
          
          message(paste("Extracting magnetometery to long format..."))
          
          # Prepare parameters for parallel processing
          params_list <- lapply(seq_len(nrow(mag)), function(i) {
            list(
              row_id = i,
              raw_values = mag$mag.magnetic.fields.raw[i],
              axes = mag$eobs.magnetometery.axes[i]
            )
          })
          
          # Process magnetometer data to long format
          processed_mag <- pblapply(params_list, process_data)
          # Combine results into a single data.table
          processed_mag <- rbindlist(processed_mag, fill = TRUE)
          
          # Add metadata back using the row_id
          processed_mag <- merge(processed_mag, mag[, .(row_id = .I, timestamp,
                                                        individual.taxon.canonical.name, tag.local.identifier, 
                                                        individual.local.identifier, eobs.magnetometery.axes, 
                                                        mag.magnetic.field.sampling.frequency.per.axis)], 
                                 by = "row_id", all.x = TRUE)
          
          # Remove redundant mag data.table
          rm(mag)
          # Identify and dynamically handle present axes
          present_axes <- intersect(c("X", "Y", "Z"), names(processed_mag))
          setcolorder(processed_mag, c(setdiff(names(processed_mag), present_axes), present_axes)) # All other columns except the present axes
          # Rename the present axes to `mag_x`, `mag_y`, and `mag_z`
          new_axis_names <- paste0("mag_", tolower(present_axes))  # Convert to mag_x, mag_y, mag_z
          setnames(processed_mag, old = present_axes, new = new_axis_names)
          
          message("Assigning burst lengths...")
          processed_mag[, burst_length := .N, by = row_id]
          
          ###############################
          message("Extrapolating within burst timestamps...")
          # Extrapolate timestamps within each burst using row_id
          processed_mag[, extrapolated_timestamp := timestamp]  # Start with the original timestamp
          
          # Add time increment based on the sampling frequency
          processed_mag[, time_increment := 1 / mag.magnetic.field.sampling.frequency.per.axis, by = row_id]
          
          # Initialize progress counter
          total_rows <- length(unique(processed_mag$row_id))
          progress_step <- floor(total_rows / 100)  # Update every 1% of progress
          message("Adjusting extrapolated timestamps...")
          
          processed_mag[, extrapolated_timestamp := {
            # Use a static counter to track progress
            if (.GRP %% progress_step == 0) {
              message(paste0("Progress: ", round(.GRP / total_rows * 100), "%"))
            }
            # Handle single sensor type dataset
            timestamp + (seq_len(.N) - 1) * time_increment
          }, by = row_id]
          
          # Drop intermediate columns if necessary
          processed_mag[, c("time_increment") := NULL]
          
          # Add a time difference column (not grouped by row_id)
          processed_mag[, time_diff := c(NA, diff(extrapolated_timestamp))]
          
          # Mark duplicates in the `extrapolated_timestamp` column
          processed_mag[, duplicate_times := duplicated(extrapolated_timestamp)]
          
          ################# Plot uninterrupted duration of each sensor type ########################
          message("Processing sampling durations and uniterupted burst lengths...")
          
          # Initialize burst_id
          processed_mag[, burst_id := {
            indices <- seq_len(.N)
            # Initialize burst_id
            burst_id <- 1
            # Create a vector to hold burst IDs for the group
            burst_ids <- integer(.N)
            burst_ids[1] <- burst_id  # Assign the first burst ID
            
            # Iterate through rows to assign burst_id based on conditions
            pbapply::pblapply(indices[-1], function(i) {
              # Assign burst_id based on conditions
              if (time_diff[i] > 0.5) { #0.5 s the cut off for continuous bursts
                burst_id <<- burst_id + 1  # Increment burst_id when conditions are met
                burst_ids[i] <<- burst_id
              } else {
                burst_ids[i] <<- burst_id  # Keep the same burst_id
              }
            })
            
            burst_ids  # Return the vector of burst IDs
            
          }]
          
        }
      }
    }
    
    ################################################################################################################################################################################################
    ################################################################################################################################################################################################
    #################################################QUATERNIONS####################################################################################################################################
    
    # Check if "orientation" exists in the sensor.type column
    if ("orientation" %in% df$sensor.type) {
      required_columns <- c(
        "timestamp", 
        "individual.taxon.canonical.name", 
        "tag.local.identifier", 
        "individual.local.identifier",
        "orientation.quaternions.sampling.frequency", 
        "orientation.quaternions.raw"
      )
      if (validate_columns(df, required_columns, "Quaternion")) {
        # Check if there are any non-NA space-separated character strings in orientation.quaternions.raw
        if (any(!is.na(df$orientation.quaternions.raw) & nzchar(df$orientation.quaternions.raw) & !is.na(df$orientation.quaternions.sampling.frequency))) {
          
          message("Starting the quaternion data processing pipeline...")
          # Filter rows for orientation data
          quat <- df[sensor.type == "orientation"]
          
          # Arrange rows by timestamp
          quat <- quat[order(timestamp)]
          rownames(quat) <- NULL
          # Make an XYZ column
          axes = "WXYZ"
          quat[, eobs.quaternion.axes := axes] 
          
          # Select and reorder
          required_columns <- c(
            "timestamp", 
            "individual.taxon.canonical.name", 
            "tag.local.identifier", 
            "individual.local.identifier",
            "eobs.quaternion.axes",
            "orientation.quaternions.sampling.frequency", 
            "orientation.quaternions.raw"
          )
          quat <- quat[, required_columns[required_columns %in% names(quat)], with = FALSE]
          setcolorder(quat, required_columns[required_columns %in% names(quat)])
          
          message(paste("Extracting quaternions to long format..."))
          
          # Function to expand raw IMU space-seperated data into long format
          process_quaternion_data <- function(params) {
            row_id <- params$row_id
            raw_values <- params$raw_values
            axes <- params$axes
            
            if (is.na(raw_values) || !is.character(raw_values)) {
              return(data.table(row_id = row_id, W = NA, X = NA, Y = NA, Z = NA))
            }
            
            # Split raw values into a numeric vector
            values <- as.numeric(unlist(strsplit(raw_values, " ")))
            
            # Determine axis order based on axes
            axis_order <- switch(
              axes,
              "WXYZ" = c("W", "X", "Y", "Z"),
            )
            
            # Repeat axis order to match the length of the raw values
            repeat_order <- rep(axis_order, length.out = length(values))
            
            # Create a data.table with row_id and expanded acceleration data
            data.table(
              row_id = row_id,
              W = if ("W" %in% repeat_order) values[repeat_order == "W"] else NA,
              X = if ("X" %in% repeat_order) values[repeat_order == "X"] else NA,
              Y = if ("Y" %in% repeat_order) values[repeat_order == "Y"] else NA,
              Z = if ("Z" %in% repeat_order) values[repeat_order == "Z"] else NA
            )
          }
          
          # Prepare parameters for parallel processing
          params_list <- lapply(seq_len(nrow(quat)), function(i) {
            list(
              row_id = i,
              raw_values = quat$orientation.quaternions.raw[i],
              axes = quat$eobs.quaternion.axes[i]
            )
          })
          
          # Process magnetometer data to long format
          processed_quat <- pblapply(params_list, process_quaternion_data)
          # Combine results into a single data.table
          processed_quat <- rbindlist(processed_quat, fill = TRUE)
          
          # Add metadata back using the row_id
          processed_quat <- merge(processed_quat, quat[, .(row_id = .I, timestamp,
                                                           individual.taxon.canonical.name, tag.local.identifier, 
                                                           individual.local.identifier, eobs.quaternion.axes, 
                                                           orientation.quaternions.sampling.frequency)], 
                                  by = "row_id", all.x = TRUE)
          
          # Remove redundant quat data.table
          rm(quat)
          # Identify and dynamically handle present axes
          present_axes <- intersect(c("W", "X", "Y", "Z"), names(processed_quat))
          setcolorder(processed_quat, c(setdiff(names(processed_quat), present_axes), present_axes)) # All other columns except the present axes
          # Rename the present axes to `mag_x`, `mag_y`, and `mag_z`
          new_axis_names <- paste0("quat_", tolower(present_axes))  # Convert to mag_x, mag_y, mag_z
          setnames(processed_quat, old = present_axes, new = new_axis_names)
          
          message("Converting quaternions...")
          # Convert quaternions
          # Step 1: Calculate the norm of the vector components
          processed_quat[, r := sqrt(quat_x^2 + quat_y^2 + quat_z^2)]
          # Step 2: Normalize quat_w
          processed_quat[, quat_w := quat_w / 32768]
          # Step 3: Calculate scale factor 's', avoiding division by zero
          processed_quat[, s := fifelse(r != 0, sqrt((1 - quat_w^2) / r), 0)]
          # Step 4: Normalize the vector components using 's'
          processed_quat[, `:=`(
            quat_x = s * quat_x,
            quat_y = s * quat_y,
            quat_z = s * quat_z
          )]
          # Step 5: Drop intermediate variables 'r' and 's'
          processed_quat[, c("r", "s") := NULL]
          
          message("Assigning burst lengths...")
          processed_quat[, burst_length := .N, by = row_id]
          
          ###############################
          
          message("Extrapolating within burst timestamps...")
          # Extrapolate timestamps within each burst using row_id
          processed_quat[, extrapolated_timestamp := timestamp]  # Start with the original timestamp
          
          # Add time increment based on the sampling frequency
          processed_quat[, time_increment := 1 / orientation.quaternions.sampling.frequency, by = row_id]
          
          # Initialize progress counter
          total_rows <- length(unique(processed_quat$row_id))
          progress_step <- floor(total_rows / 100)  # Update every 1% of progress
          message("Adjusting extrapolated timestamps...")
          
          processed_quat[, extrapolated_timestamp := {
            # Use a static counter to track progress
            if (.GRP %% progress_step == 0) {
              message(paste0("Progress: ", round(.GRP / total_rows * 100), "%"))
            }
            # Handle single sensor type dataset
            timestamp + (seq_len(.N) - 1) * time_increment
          }, by = row_id]
          
          # Drop intermediate columns if necessary
          processed_quat[, c("time_increment") := NULL]
          
          # Add a time difference column (not grouped by row_id)
          processed_quat[, time_diff := c(NA, diff(extrapolated_timestamp))]
          
          # Mark duplicates in the `extrapolated_timestamp` column
          processed_quat[, duplicate_times := duplicated(extrapolated_timestamp)]
          
          ################# Plot uninterrupted duration of each sensor type ########################
          message("Processing sampling durations and uniterupted burst lengths...")
          
          # Initialize burst_id
          processed_quat[, burst_id := {
            indices <- seq_len(.N)
            # Initialize burst_id
            burst_id <- 1
            # Create a vector to hold burst IDs for the group
            burst_ids <- integer(.N)
            burst_ids[1] <- burst_id  # Assign the first burst ID
            
            # Iterate through rows to assign burst_id based on conditions
            pbapply::pblapply(indices[-1], function(i) {
              # Assign burst_id based on conditions
              if (time_diff[i] > 0.5) { #0.5 s the cut off for continuous bursts
                burst_id <<- burst_id + 1  # Increment burst_id when conditions are met
                burst_ids[i] <<- burst_id
              } else {
                burst_ids[i] <<- burst_id  # Keep the same burst_id
              }
            })
            
            burst_ids  # Return the vector of burst IDs
            
          }]
          
        }
      }
    }
    
    # Remove df
    rm(df, required_columns) 
    
    # Final processing logic to check which data.tables exist and handle accordingly
    result <- NULL
    
    # Check which processed data.tables exist
    processed_acc_exists <- exists("processed_acc")
    processed_mag_exists <- exists("processed_mag")
    processed_quat_exists <- exists("processed_quat")
    
    if (processed_acc_exists && processed_mag_exists && processed_quat_exists) {
      # Case 1: All three data.tables exist
      if (nrow(processed_mag) == nrow(processed_quat) &&
          all(processed_mag$extrapolated_timestamp == processed_quat$extrapolated_timestamp)) {
        # processed_mag and processed_quat have matching row numbers and timestamps
        message("Combining Magnetometer and Quaternion data. Acceleration & Orientation data saved in list...")
        combined_data <- cbind(
          processed_mag,
          processed_quat[, .(quat_w = quat_w, quat_x = quat_x, quat_y = quat_y, quat_z = quat_z)]
        )
        result <- list(
          `Acceleration Data` = processed_acc,
          `Orientation Data` = combined_data
        )
      } else {
        # processed_mag and processed_quat do not match
        message("Acceleration, Magnetometer and Quaternion data saved in list...")
        result <- list(
          `Acceleration Data` = processed_acc,
          `Magnetometery Data` = processed_mag,
          `Quaternion Data` = processed_quat
        )
      }
    } else if (processed_acc_exists && processed_mag_exists) {
      # Case 3a: Only Acceleration and Magnetometer data exist
      message("Acceleration and Magnetometer data saved in list...")
      result <- list(
        `Acceleration Data` = processed_acc,
        `Magnetometery Data` = processed_mag
      )
    } else if (processed_acc_exists && processed_quat_exists) {
      # Case 3b: Only Acceleration and Quaternion data exist
      message("Acceleration and Quaternion data saved in list...")
      result <- list(
        `Acceleration Data` = processed_acc,
        `Quaternion Data` = processed_quat
      )
    } else if (processed_mag_exists && processed_quat_exists) {
      # Case 4: Only Magnetometer and Quaternion data exist
      if (nrow(processed_mag) == nrow(processed_quat) &&
          all(processed_mag$extrapolated_timestamp == processed_quat$extrapolated_timestamp)) {
        # processed_mag and processed_quat have matching row numbers and timestamps
        message("No Acceleration data. Combining Magnetometer and Quaternion data into a single data frame...")
        combined_data <- cbind(
          processed_mag,
          processed_quat[, .(quat_w = quat_w, quat_x = quat_x, quat_y = quat_y, quat_z = quat_z)]
        )
        result <- as.data.frame(`Orientation Data` = combined_data)
      } else {
        message("Magnetometer and Quaternion data saved in list...")
        # processed_mag and processed_quat do not match
        result <- list(
          `Magnetometery Data` = processed_mag,
          `Quaternion Data` = processed_quat
        )
      }
    } else if (processed_acc_exists) {
      message("Only Acceleration data present. Data saved as a single data frame...")
      # Case 5: Only Acceleration data exists
      result <- as.data.frame(processed_acc)
    } else {
      # No valid processed data exists
      stop("No processed data is available to return.")
    }
    
    # Return the result
    return(result)
    gc()
    
  })
  
}

#################################### END OF FUNCTION ########################################
###########################################################################################################################################################################################################   
########################################################################################################################################################################################################### 

# E.g., 

#Animal.1 = Eobs_Data_Reader(df = df, 
                        # rolling_mean_width = 40, 
                        # standardised_freq_rate = NULL, 
                        # start_timestamp = NULL,
                        # end_timestamp = NULL,
                        # plot = TRUE)
