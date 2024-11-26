setwd("/Users/zinomeyer/Documents/06-coding/image-analysis-arctic-icecaps/")

install.packages("imager")
library(imager)

icecap_1979 <- load.image("data/sea-ice-concentration-map-1979-aug.jpeg")
icecap_1989 <- load.image("data/sea-ice-concentration-map-1989-aug.jpeg")
icecap_1999 <- load.image("data/sea-ice-concentration-map-1999-aug.jpeg")
icecap_2009 <- load.image("data/sea-ice-concentration-map-2009-aug.jpeg")
icecap_2019 <- load.image("data/sea-ice-concentration-map-2019-aug.jpeg")
# dim(icecap_1979)

icecap_1979_gray <- grayscale(icecap_1979)
icecap_1989_gray <- grayscale(icecap_1989)
icecap_1999_gray <- grayscale(icecap_1999)
icecap_2009_gray <- grayscale(icecap_2009)
icecap_2019_gray <- grayscale(icecap_2019)
plot(icecap_1989_gray)

# SVD analysis
icecap_1979_svd <- svd(icecap_1979_gray)
icecap_1989_svd <- svd(icecap_1989_gray)
icecap_1999_svd <- svd(icecap_1999_gray)
icecap_2009_svd <- svd(icecap_2009_gray)
icecap_2019_svd <- svd(icecap_2019_gray)

# Print singular values, EOFs and PCs for 1979
icecap_1979_svd$d[1:5]
icecap_1979_svd$u[1:5, 1:5]
icecap_1979_svd$v[1:5, 1:5]

# Print singular values, EOFs and PCs for 1989
icecap_1989_svd$d[1:5]
icecap_1989_svd$u[1:5, 1:5]
icecap_1989_svd$v[1:5, 1:5]

# Print singular values, EOFs and PCs for 1999
icecap_1999_svd$d[1:5]
icecap_1999_svd$u[1:5, 1:5]
icecap_1999_svd$v[1:5, 1:5]

# Print singular values, EOFs and PCs for 2009
icecap_2009_svd$d[1:5]
icecap_2009_svd$u[1:5, 1:5]
icecap_2009_svd$v[1:5, 1:5]

# Print singular values, EOFs and PCs for 2019
icecap_2019_svd$d[1:5]
icecap_2019_svd$u[1:5, 1:5]
icecap_2019_svd$v[1:5, 1:5]


# Print all the singular values of all years in a table
t <- 1:10
all_singular_values <- data.frame(
    icecap_1979_svd$d[t],
    icecap_1989_svd$d[t],
    icecap_1999_svd$d[t],
    icecap_2009_svd$d[t],
    icecap_2019_svd$d[t]
)
colnames(all_singular_values) <- c("1979", "1989", "1999", "2009", "2019")
all_singular_values



# Calculate means of the singular values for each year
mean_singular_values <- data.frame(
    Year = c("1979", "1989", "1999", "2009", "2019"),
    Mean_Singular_Value = c(
        mean(icecap_1979_svd$d[t]),
        mean(icecap_1989_svd$d[t]),
        mean(icecap_1999_svd$d[t]),
        mean(icecap_2009_svd$d[t]),
        mean(icecap_2019_svd$d[t])
    )
)

# Print the means of the singular values
print(mean_singular_values)

# Create a bar plot of the mean singular values
barplot(mean_singular_values$Mean_Singular_Value,
    names.arg = mean_singular_values$Year,
    col = "skyblue",
    main = "Mean Singular Values of Icecap Images by Year",
    xlab = "Year",
    ylab = "Mean Singular Value",
    ylim = c(0, max(mean_singular_values$Mean_Singular_Value) * 1.2)
)


# Plot singular values
variance_percent_D <- 100 * (icecap_1979_svd$d^2) / sum(icecap_1979_svd$d^2)
cum_percent_D <- cumsum(variance_percent_D)
modeK <- 1:length(icecap_1979_svd$d)
K <- 30

plot(modeK[1:K], variance_percent_D[1:K],
    type = "o", col = "blue",
    xlab = "Mode number", pch = 16,
    ylab = "Percentage of mode variance",
    main = "Scree Plot of first 30 eigenvalues of B/W Icecap Image 1979"
)

dev.off()


# Scree plot of first 30 modes of A, incl. cumulative variance
par(mar = c(4, 4, 2, 4), mgp = c(2.2, 0.7, 0))
plot(1:K,
    variance_percent_D[1:K],
    ylim = c(0, 100),
    type = "o",
    ylab = "Percentage of Variance [%]",
    xlab = "EOF Mode Number",
    cex.lab = 1.2, cex.axis = 1.1, lwd = 2,
    main = "Scree Plot of first 30 eigenvalues of B/W Sunset Photo"
)
legend(3, 30,
    col = c("black"), lty = 1, lwd = 2.0,
    legend = c("Percentange Variance"), bty = "n",
    text.font = 2, cex = 1.0, text.col = "black"
)

par(new = TRUE)
plot(1:K, cum_percent_D[1:K],
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






# Calculate variance and cumulative variance for each year
calculate_variance <- function(svd_result) {
    variance_percent <- 100 * (svd_result$d^2) / sum(svd_result$d^2)
    cum_percent <- cumsum(variance_percent)
    return(list(variance_percent = variance_percent, cum_percent = cum_percent))
}

variance_1979 <- calculate_variance(icecap_1979_svd)
variance_1989 <- calculate_variance(icecap_1989_svd)
variance_1999 <- calculate_variance(icecap_1999_svd)
variance_2009 <- calculate_variance(icecap_2009_svd)
variance_2019 <- calculate_variance(icecap_2019_svd)

# Plot scree plot with variance and cumulative variance for each year
K <- 30
par(mar = c(4, 4, 2, 4), mgp = c(2.2, 0.7, 0))
plot(1:K, variance_1979$variance_percent[1:K],
    ylim = c(0, 100), type = "o", ylab = "Percentage of Variance [%]",
    xlab = "EOF Mode Number", cex.lab = 1.2, cex.axis = 1.1, lwd = 2,
    main = "Scree Plot of first 30 eigenvalues of B/W Icecap Images"
)
lines(1:K, variance_1989$variance_percent[1:K], type = "o", col = "red", lwd = 2)
lines(1:K, variance_1999$variance_percent[1:K], type = "o", col = "green", lwd = 2)
lines(1:K, variance_2009$variance_percent[1:K], type = "o", col = "purple", lwd = 2)
lines(1:K, variance_2019$variance_percent[1:K], type = "o", col = "orange", lwd = 2)

legend("topright",
    legend = c("1979", "1989", "1999", "2009", "2019"),
    col = c("black", "red", "green", "purple", "orange"),
    lty = 1, lwd = 2, bty = "n"
)

par(new = TRUE)
plot(1:K, variance_1979$cum_percent[1:K],
    ylim = c(90, 100), type = "o", col = "blue", lwd = 2, axes = FALSE,
    xlab = "", ylab = ""
)
lines(1:K, variance_1989$cum_percent[1:K], type = "o", col = "red", lwd = 2)
lines(1:K, variance_1999$cum_percent[1:K], type = "o", col = "green", lwd = 2)
lines(1:K, variance_2009$cum_percent[1:K], type = "o", col = "purple", lwd = 2)
lines(1:K, variance_2019$cum_percent[1:K], type = "o", col = "orange", lwd = 2)

legend("bottomright",
    legend = c("Cumulative 1979", "Cumulative 1989", "Cumulative 1999", "Cumulative 2009", "Cumulative 2019"),
    col = c("blue", "red", "green", "purple", "orange"),
    lty = 1, lwd = 2, bty = "n"
)

axis(4, col = "blue", col.axis = "blue", mgp = c(3, 0.7, 0))
mtext("Cumulative Variance [%]", col = "blue", cex = 1.2, side = 4, line = 2)







# Load necessary libraries
install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")
library(readxl)
library(dplyr)
library(tidyr)

setwd("~/Documents/04-projects/image-analysis-arctic-icecaps")


# Define the path to the Excel file
file_path <- "data/N_Sea_Ice_Index_Regional_Monthly_Data_G02135_v3.0.xlsx"

# Get all sheet names
sheet_names <- excel_sheets(file_path)

# Separate sheets into area and extent
area_sheets <- sheet_names[grepl("area", sheet_names, ignore.case = TRUE)]
extent_sheets <- sheet_names[grepl("extent", sheet_names, ignore.case = TRUE)]

months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")

process_sheets <- function(sheets) {
    all_regions <- list()

    for (sheet in sheets) {
        xls_data <- read_excel(file_path, sheet = sheet)
        region_data_df <- data.frame(xls_data)

        # Get years
        years <- region_data_df[3:49, 1]

        # Remove inimportant rows and columns
        region_area <- data.frame(region_data_df[3:49, seq(2, ncol(region_data_df), by = 2)])

        # Cell formatting (numeric, no NA, 2 decimal points)
        region_area[is.na(region_area)] <- 0
        region_area <- as.data.frame(lapply(region_area, as.numeric))
        region_area <- round(region_area, 2)

        # Set row and column names
        rownames(region_area) <- years
        colnames(region_area) <- months

        # Flatten the data into one row
        region_space_time_df <- reshape(cbind(region_area, id = 1, time = rownames(region_area)),
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

# Process and merge all area sheets
area_space_time_df <- process_sheets(area_sheets)

# Process and merge all extent sheets
extent_space_time_df <- process_sheets(extent_sheets)

# # Round all cells to 2 decimal points
# area_space_time_df <- round(area_space_time_df, 2)
# extent_space_time_df <- round(extent_space_time_df, 2)

# Display an excerpt of the final data frames
print(dim(area_space_time_df))
print(area_space_time_df[, c("January-1978", "August-2023")])

print(dim(extent_space_time_df))
print(extent_space_time_df[, c("August-1980", "August-2023")])
