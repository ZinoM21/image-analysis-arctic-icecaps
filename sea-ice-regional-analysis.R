# Load necessary libraries
install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggOceanMaps")
install.packages("ggspatial")

library(readxl)
library(dplyr)
library(tidyr)
library(ggOceanMaps)
library(ggspatial)

# Define the path to the Excel file
file_path <- "data/N_Sea_Ice_Index_Regional_Monthly_Data_G02135_v3.0.xlsx"

# Get all sheet names
sheet_names <- excel_sheets(file_path)

# Separate sheets into area and extent
area_sheets <- sheet_names[grepl("area", sheet_names, ignore.case = TRUE)]
extent_sheets <- sheet_names[grepl("extent", sheet_names, ignore.case = TRUE)]

months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")

process_sheets <- function(sheets, avg_rows) {
    all_regions <- list()

    for (sheet in sheets) {
        xls_data <- read_excel(file_path, sheet = sheet)
        region_data_df <- data.frame(xls_data)

        # Get years
        years <- region_data_df[3:49, 1]

        # Remove inimportant rows and columns
        region_area <- data.frame(region_data_df[3:49, seq(2, ncol(region_data_df), by = 2)])

        # Cell formatting (numeric, no NA, 2 decimal points)
        region_area <- as.data.frame(lapply(region_area, as.numeric))
        region_area <- round(region_area, 2)

        # Set row and column names
        rownames(region_area) <- years
        colnames(region_area) <- months

        # Replace values with values minus each colmeans
        means <- colMeans(region_area[avg_rows, ], na.rm = TRUE)
        region_area_anomalies <- sweep(region_area, 2, means, "-")

        # Remove NA values
        region_area_anomalies[is.na(region_area_anomalies)] <- 0

        # Flatten the data into one row
        region_space_time_df <- reshape(cbind(region_area_anomalies, id = 1, time = rownames(region_area_anomalies)),
            direction = "wide", sep = "-", new.row.names = 0
        )[-1]

        # Format region name
        region_name <- sub("-(Area|Extent).*", "", sheet)
        region_name <- gsub("-", " ", region_name)

        # Add region row to list
        all_regions[[region_name]] <- region_space_time_df
    }

    all_regions_dataframe <- do.call(rbind, all_regions)
    return(all_regions_dataframe)
}

# Set median timeframe in years, starts at 1979 (1 = 1979)
rows_for_mean_calc <- 1:20

# Process and merge all area sheets
area_space_time_df <- process_sheets(area_sheets, rows_for_mean_calc)

# Process and merge all extent sheets
extent_space_time_df <- process_sheets(extent_sheets, rows_for_mean_calc)

# Display an excerpt of the final data frames
print(dim(area_space_time_df))
print(area_space_time_df[, c("January-1980", "August-2023")])

print(dim(extent_space_time_df))
print(extent_space_time_df[, c("August-1980", "August-2023")])

# Export the data frames to CSV
# write.csv(area_space_time_df, "data/sea-ice-regional-area-anomalies.csv")
# write.csv(extent_space_time_df, "data/sea-ice-regional-extent-anomalies.csv")


# Visualize data ------------
# Heat map
generate_heat_map <- function(data, title, month = "August") {
    month_data <- data[, grepl(month, colnames(data))]
    month_data <- month_data[, seq(3, ncol(month_data), by = 10)]

    region_names <- rownames(month_data)
    time_stamps <- gsub(".*-(\\d{4})", "\\1", colnames(month_data))

    par(mar = c(4, 10, 4, 4) + 0.1)
    image(t(month_data),
        col = heat.colors(64),
        axes = FALSE
    )
    axis(1, at = seq(0, 1, length.out = 5), labels = time_stamps)
    axis(2, at = seq(0, 1, length.out = 14), labels = region_names, las = 2)
    title(main = title)
    box()
}

generate_heat_map(
    area_space_time_df,
    "Regional sea ice area anomalies in September over time [sqkm]",
    "September"
)

generate_heat_map(
    extent_space_time_df,
    "Regional sea ice extent anomalies in September over time [sqkm]",
    "September"
)


# SVD analysis ---------------
svd_area <- svd(area_space_time_df)
U_area <- svd_area$u
D_area <- svd_area$d
V_area <- svd_area$v
EOF1_area <- U_area[1, ]

svd_extent <- svd(extent_space_time_df)
U_extent <- svd_extent$u
D_extent <- svd_extent$d
V_extent <- svd_extent$v
EOF1_extent <- U_extent[1, ]


# Display the first EOF for each matrix
print(U_area[, 1])
print(U_extent[, 1])

SVDd <- svd_area$d

# Scree plot ---
# Set needed values
lam <- SVDd^2
K <- 20
lamK <- lam[1:K]

percentD <- 100 * (lam) / sum(lam)
cumpercentD <- cumsum(percentD)
modeK <- 1:length(SVDd)

# Scree Plot
plot(modeK[1:K], percentD[1:K],
    type = "o", col = "blue",
    xlab = "Mode number", pch = 16,
    ylab = "Percentage of mode variance",
    main = "Scree Plot of SVD"
)

# Scree Plot with Cumulative Percentage
par(mar = c(4, 4, 2, 4), mgp = c(2.2, 0.7, 0))
plot(1:K, 100 * lamK / sum(lam),
    ylim = c(0, 100), type = "o",
    ylab = "Percentage of Variance [%]",
    xlab = "EOF Mode Number",
    cex.lab = 1.2, cex.axis = 1.1, lwd = 2,
    main = "Scree Plot of the First 20 Eigenvalues"
)
legend(3, 30,
    col = c("black"), lty = 1, lwd = 2.0,
    legend = c("Percentange Variance"), bty = "n",
    text.font = 2, cex = 1.0, text.col = "black"
)
par(new = TRUE)
plot(1:K, cumsum(100 * lamK / sum(lam)),
    ylim = c(90, 100), type = "o",
    col = "blue", lwd = 2, axes = FALSE,
    xlab = "", ylab = ""
)
legend(3, 94.5,
    col = c("blue"), lty = 1, lwd = 2.0,
    legend = c("Cumulative Percentage Variance"), bty = "n",
    text.font = 2, cex = 1.0, text.col = "blue"
)
axis(4, col = "blue", col.axis = "blue", mgp = c(3, 0.7, 0))
mtext("Cumulative Variance [%]",
    col = "blue",
    cex = 1.2, side = 4, line = 2
)


# -------- Plot EOF mode 1 on Map --------
# Run those lines at top of file:
# install.packages("ggOceanMaps")
# install.packages("ggspatial")
# library(ggOceanMaps)
# library(ggspatial)

# Coords (by Google) & Order for all regions
# baffin = 72.809094, -66.175709
# barents = 74.576217, 37.519449
# beaufort = 72.961356, -145.563116
# bering = 56.855056, -177.762477
# canadian = 74.449727, -102.979623
# central arctic = 90, 0
# chukchi = 69.463991, -171.928849
# east siberian = 73.817176, 157.271205
# greenland = 77.850046, -5.159874
# hudson = 59.710451, -85.768872
# kara = 75.151568, 73.153283
# laptec = 76.032029, 126.462445
# okhotsk = 52.978659, 149.296564
# st lawrence = 48.552296, -61.379056


# Combine data for easy access
locations <- data.frame(
    name = rownames(area_space_time_df),
    lat = c(
        72.809094, 74.576217, 72.961356,
        56.855056, 74.449727, 90,
        69.463991, 73.817176, 77.850046,
        59.710451, 75.151568, 76.032029,
        52.978659, 48.552296
    ),
    lon = c(
        -66.175709, 37.519449, -145.563116,
        -177.762477, -102.979623, 0,
        -171.928849, 157.271205, -5.159874,
        -85.768872, 73.153283, 126.462445,
        149.296564, -61.379056
    ),
    EOF_value = EOF1_area,
    hjust = c(0.5, 0.5, 0.8, 0.5, 0.5, 0.5, 0.5, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.3),
    vjust = c(-0.5, -0.5, -0.5, -0.8, -0.7, -0.6, -0.5, -0.5, 1.7, -0.5, -0.5, -0.5, -0.5, -0.5)
)

# Set color dependant on value
locations$cols <- ifelse(locations$EOF_value > 0, "blue", "red")

# Plot Mode 1
basemap(
    limits = 45,
    shapefiles = "Arctic",
    bathymetry = TRUE, # optional
    bathy.style = "raster_binned_blues", # optional
    land.col = "#eeeac4", # optional
    legends = FALSE
) +
    geom_point(
        data = transform_coord(locations),
        aes(x = lon, y = lat, size = abs(locations$EOF_value)),
        color = locations$cols,
        show.legend = FALSE
    ) +
    geom_label(
        data = transform_coord(locations),
        aes(x = lon, y = lat),
        label = paste(locations$name, round(locations$EOF_value, 3)),
        hjust = locations$hjust,
        vjust = locations$vjust,
        size = 3,
        color = locations$col
    )

dev.off()
