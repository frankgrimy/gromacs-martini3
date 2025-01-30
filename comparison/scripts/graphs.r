# Check if the Peptides package is installed, if not install it
if (!requireNamespace("Peptides", quietly = TRUE)) {
  install.packages("Peptides")
}
library(Peptides)
library(ggplot2)
options(digits = 10)

REPOROOT <- system("git rev-parse --show-toplevel", intern = TRUE)
BASEDIR <- (paste(REPOROOT, "/comparison/workdir/base", sep = ""))
RESCALDIR <- (paste(REPOROOT, "/comparison/workdir/rescaled", sep = ""))
GRAPHDIR <- (paste(REPOROOT, "/comparison/graphs", sep = ""))

### POTENTIAL ENERGY PLOT (Minimization) ###
potbase <- readXVG(paste(BASEDIR, "/potential.xvg", sep = ""))
potrescaled <- readXVG(paste(RESCALDIR, "/potential.xvg", sep = ""))

# Transform the data frame column data to numeric
potbase$Time <- as.numeric(potbase$Time)
potbase$Potential <- as.numeric(potbase$Potential)
potrescaled$Time <- as.numeric(potrescaled$Time)
potrescaled$Potential <- as.numeric(potrescaled$Potential)

# Transform time from ps to ns
potbase$Time <- potbase$Time / 1000
potrescaled$Time <- potrescaled$Time / 1000

# Plot the potential energy

pot_plot <- ggplot(potbase, aes(x = Time)) +
  geom_line(aes(y = Potential, color = "Base"), linewidth = 0.5) +
  geom_line(data = potrescaled, aes(y = Potential, color = "Rescaled"), linewidth = 0.5) +
  scale_color_manual(values = c("Base" = "blue", "Rescaled" = "#FF0000C0"),
                    labels = c("Base", "Reescalado"),
                    name = "Series") +
  scale_y_continuous(limits = c(-1.5e6, -7e5), breaks = seq(-1.5e6, -7e5, 1.25e5)) +
  labs(title = "Energía Potencial durante Minimización",
      x = "Tiempo (ns)",
      y = "Energía Potencial (kJ/mol)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(paste(GRAPHDIR, "/minim_pot.svg", sep = ""), pot_plot, width = 8, height = 6, dpi = 300)
#print(pot_plot)

### DENSITY PLOT ###
density <- readXVG(paste(BASEDIR, "/density.xvg", sep = ""))

# Apend the rescaled density to the density data frame
density$Rescaled <- readXVG(paste(RESCALDIR, "/density.xvg", sep = ""))$Density

# Transform the data frame column data to numeric
density$Time <- as.numeric(density$Time)
density$Density <- as.numeric(density$Density)
density$Rescaled <- as.numeric(density$Rescaled)
# Transform time from ps to ns
density$Time <- density$Time / 1000

# Calculate mean and sd for the base and rescaled densities
mean_base <- signif(mean(density$Density), 6)
mean_rescaled <- signif(mean(density$Rescaled), 6)
sd_base <- signif(sd(density$Density), 3)
sd_rescaled <- signif(sd(density$Rescaled), 3)

# Mean and sd for the first 0.5 ns
mean_base_05ns <- signif(mean(density$Density[density$Time <= 0.5]), 6)
mean_rescaled_05ns <- signif(mean(density$Rescaled[density$Time <= 0.5]), 6)
sd_base_05ns <- signif(sd(density$Density[density$Time <= 0.5]), 3)
sd_rescaled_05ns <- signif(sd(density$Rescaled[density$Time <= 0.5]), 3)

# Print the mean and sd for the base and rescaled densities
cat("Densidades (kg/m³):\n")
cat("Base: ", mean_base, " ± ", sd_base, "\n")
cat("Reescalado: ", mean_rescaled, " ± ", sd_rescaled, "\n")
cat("Base (0.5 ns): ", mean_base_05ns, " ± ", sd_base_05ns, "\n")
cat("Reescalado (0.5 ns): ", mean_rescaled_05ns, " ± ", sd_rescaled_05ns, "\n")
cat("\n")

# Plot both densities
dens_plot <- ggplot(density, aes(x = Time)) +
  geom_line(aes(y = Density, color = "Density"), linewidth = 0.5) +
  geom_line(aes(y = Rescaled, color = "Rescaled"), linewidth = 0.5) +
  scale_color_manual(values = c("Density" = "blue", "Rescaled" = "#FF0000C0"),
                     labels = c("Base", "Reescalado"),
                     name = "Series") +
  scale_y_continuous(limits = c(950, 1050), breaks = seq(950, 1050, 25)) +
  labs(title = "Densidad durante Equilibración NPT",
       x = "Tiempo (ns)",
       y = "Densidad (kg/m³)") +
  theme_minimal() +
  theme(legend.position = "right")

dens_plot_05ns <- dens_plot + scale_x_continuous(limits = c(0, 0.5)) +
  labs(title = "Densidad durante Equilibración NPT (0.5 ns)")

ggsave(paste(GRAPHDIR, "/density.svg", sep = ""), dens_plot, width = 8, height = 6, dpi = 300)
ggsave(paste(GRAPHDIR, "/density_05ns.svg", sep = ""), dens_plot_05ns, width = 8, height = 6, dpi = 300)
# print(dens_plot)
# print(dens_plot_05ns)

### PRESSURE PLOT ###
pressure <- readXVG(paste(BASEDIR, "/pressure.xvg", sep = ""))
pressure$Rescaled <- readXVG(paste(RESCALDIR, "/pressure.xvg", sep = ""))$Pressure

# Transform the data frame column data to numeric
pressure$Time <- as.numeric(pressure$Time)
pressure$Pressure <- as.numeric(pressure$Pressure)
pressure$Rescaled <- as.numeric(pressure$Rescaled)

# Transform time from ps to ns
pressure$Time <- pressure$Time / 1000

# Calculate mean and sd for the base and rescaled pressures
mean_base <- signif(mean(pressure$Pressure), 3)
mean_rescaled <- signif(mean(pressure$Rescaled), 3)
sd_base <- signif(sd(pressure$Pressure), 3)
sd_rescaled <- signif(sd(pressure$Rescaled), 3)

# Mean and sd for the first 0.5 ns
mean_base_05ns <- signif(mean(pressure$Pressure[pressure$Time <= 0.5]), 3)
mean_rescaled_05ns <- signif(mean(pressure$Rescaled[pressure$Time <= 0.5]), 3)
sd_base_05ns <- signif(sd(pressure$Pressure[pressure$Time <= 0.5]), 3)
sd_rescaled_05ns <- signif(sd(pressure$Rescaled[pressure$Time <= 0.5]), 3)

# Print the mean and sd for the base and rescaled pressures
cat("Presiones (bar):\n")
cat("Base: ", mean_base, " ± ", sd_base, "\n")
cat("Reescalado: ", mean_rescaled, " ± ", sd_rescaled, "\n")
cat("Base (0.5 ns): ", mean_base_05ns, " ± ", sd_base_05ns, "\n")
cat("Reescalado (0.5 ns): ", mean_rescaled_05ns, " ± ", sd_rescaled_05ns, "\n")
cat("\n")

# Plot both pressures
pres_plot <- ggplot(pressure, aes(x = Time)) +
  geom_line(aes(y = Pressure, color = "Pressure"), linewidth = 0.5) +
  geom_line(aes(y = Rescaled, color = "Rescaled"), linewidth = 0.5) +
  scale_color_manual(values = c("Pressure" = "blue", "Rescaled" = "#FF0000C0"),
                     labels = c("Base", "Reescalado"),
                     name = "Series") +
  labs(title = "Presión durante Equilibración NPT",
       x = "Tiempo (ns)",
       y = "Presión (bar)") +
  scale_y_continuous(limits = c(-2000, 500), breaks = seq(-2000, 500, 500)) +
  theme_minimal() +
  theme(legend.position = "right")

pres_plot_05ns <- pres_plot + scale_x_continuous(limits = c(0, 0.5)) +
  labs(title = "Presión durante Equilibración NPT (0.5 ns)")

ggsave(paste(GRAPHDIR, "/pressure.svg", sep = ""), pres_plot, width = 8, height = 6, dpi = 300)
ggsave(paste(GRAPHDIR, "/pressure_05ns.svg", sep = ""), pres_plot_05ns, width = 8, height = 6, dpi = 300)

# print(pres_plot)
# print(pres_plot_05ns)

### GYRATION RADIUS PLOT ###
gyrbase <- readXVG(paste(BASEDIR, "/gyrate.xvg", sep = ""))
gyrscaled <- readXVG(paste(RESCALDIR, "/gyrate.xvg", sep = ""))

# Transform the data frame column data to numeric
gyrbase$Time <- as.numeric(gyrbase$Time)
gyrbase$Rg <- as.numeric(gyrbase$Rg)
gyrbase$'Rg/sX/N' <- as.numeric(gyrbase$'Rg/sX/N')
gyrbase$'Rg/sY/N' <- as.numeric(gyrbase$'Rg/sY/N')
gyrbase$'Rg/sZ/N' <- as.numeric(gyrbase$'Rg/sZ/N')

gyrscaled$Time <- as.numeric(gyrscaled$Time)
gyrscaled$Rg <- as.numeric(gyrscaled$Rg)
gyrscaled$'Rg/sX/N' <- as.numeric(gyrscaled$'Rg/sX/N')
gyrscaled$'Rg/sY/N' <- as.numeric(gyrscaled$'Rg/sY/N')
gyrscaled$'Rg/sZ/N' <- as.numeric(gyrscaled$'Rg/sZ/N')

# Transform time from ps to ns
gyrbase$Time <- gyrbase$Time / 1000
gyrscaled$Time <- gyrscaled$Time / 1000

# Calculate mean and sd for the base and rescaled gyration radius
mean_base <- signif(mean(gyrbase$Rg), 3)
mean_rescaled <- signif(mean(gyrscaled$Rg), 3)
sd_base <- signif(sd(gyrbase$Rg), 3)
sd_rescaled <- signif(sd(gyrscaled$Rg), 3)

# Print the mean and sd for the base and rescaled gyration radius
cat("Radio de Giro (nm):\n")
cat("Base: ", mean_base, " ± ", sd_base, "\n")
cat("Reescalado: ", mean_rescaled, " ± ", sd_rescaled, "\n")
cat("\n")

# Calculate mean and sd for the base and rescaled gyration radius in X
mean_base_x <- signif(mean(gyrbase$'Rg/sX/N'), 3)
mean_rescaled_x <- signif(mean(gyrscaled$'Rg/sX/N'), 3)
sd_base_x <- signif(sd(gyrbase$'Rg/sX/N'), 3)
sd_rescaled_x <- signif(sd(gyrscaled$'Rg/sX/N'), 3)

# Print the mean and sd for the base and rescaled gyration radius in X
cat("Radio de Giro (eje X) (nm):\n")
cat("Base: ", mean_base_x, " ± ", sd_base_x, "\n")
cat("Reescalado: ", mean_rescaled_x, " ± ", sd_rescaled_x, "\n")
cat("\n")

# Calculate mean and sd for the base and rescaled gyration radius in Y
mean_base_y <- signif(mean(gyrbase$'Rg/sY/N'), 3)
mean_rescaled_y <- signif(mean(gyrscaled$'Rg/sY/N'), 3)
sd_base_y <- signif(sd(gyrbase$'Rg/sY/N'), 3)
sd_rescaled_y <- signif(sd(gyrscaled$'Rg/sY/N'), 3)

# Print the mean and sd for the base and rescaled gyration radius in Y
cat("Radio de Giro (eje Y) (nm):\n")
cat("Base: ", mean_base_y, " ± ", sd_base_y, "\n")
cat("Reescalado: ", mean_rescaled_y, " ± ", sd_rescaled_y, "\n")
cat("\n")

# Calculate mean and sd for the base and rescaled gyration radius in Z
mean_base_z <- signif(mean(gyrbase$'Rg/sZ/N'), 3)
mean_rescaled_z <- signif(mean(gyrscaled$'Rg/sZ/N'), 3)
sd_base_z <- signif(sd(gyrbase$'Rg/sZ/N'), 3)
sd_rescaled_z <- signif(sd(gyrscaled$'Rg/sZ/N'), 3)

# Print the mean and sd for the base and rescaled gyration radius in Z
cat("Radio de Giro (eje Z) (nm):\n")
cat("Base: ", mean_base_z, " ± ", sd_base_z, "\n")
cat("Reescalado: ", mean_rescaled_z, " ± ", sd_rescaled_z, "\n")
cat("\n")

# Plot general gyration radius
gyr_plot <- ggplot(gyrbase, aes(x = Time)) +
  geom_line(aes(y = Rg, color = "Base"), linewidth = 0.5) +
  geom_line(data = gyrscaled, aes(y = Rg, color = "Rescaled"), linewidth = 0.5) +
  scale_color_manual(values = c("Base" = "blue", "Rescaled" = "#FF0000C0"),
                     labels = c("Base", "Reescalado"),
                     name = "Series") +
  scale_y_continuous(limits = c(1, 6), breaks = seq(1, 6, 1)) +
  labs(title = "Radio de Giro durante Producción",
       x = "Tiempo (ns)",
       y = "Radio de Giro (nm)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(paste(GRAPHDIR, "/gyr.svg", sep = ""), gyr_plot, width = 8, height = 6, dpi = 300)
#print(gyr_plot)

# Plot X gyration radius
gyrx_plot <- ggplot(gyrbase, aes(x = Time)) +
  geom_line(aes(y = `Rg/sX/N`, color = "Base"), linewidth = 0.5) +
  geom_line(data = gyrscaled, aes(y = `Rg/sX/N`, color = "Rescaled"), linewidth = 0.5) +
  scale_color_manual(values = c("Base" = "blue", "Rescaled" = "#FF0000C0"),
                     labels = c("Base", "Reescalado"),
                     name = "Series") +
  scale_y_continuous(limits = c(1, 6), breaks = seq(1, 6, 1)) +
  labs(title = "Radio de Giro (eje X) durante Producción",
       x = "Tiempo (ns)",
       y = "Radio de Giro en X (nm)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(paste(GRAPHDIR, "/gyr_x.svg", sep = ""), gyrx_plot, width = 8, height = 6, dpi = 300)
# print(gyrx_plot)

# Plot Y gyration radius
gyry_plot <- ggplot(gyrbase, aes(x = Time)) +
  geom_line(aes(y = `Rg/sY/N`, color = "Base"), linewidth = 0.5) +
  geom_line(data = gyrscaled, aes(y = `Rg/sY/N`, color = "Rescaled"), linewidth = 0.5) +
  scale_color_manual(values = c("Base" = "blue", "Rescaled" = "#FF0000C0"),
                     labels = c("Base", "Reescalado"),
                     name = "Series") +
  scale_y_continuous(limits = c(1, 6), breaks = seq(1, 6, 1)) +
  labs(title = "Radio de Giro (eje Y) durante Producción",
       x = "Tiempo (ns)",
       y = "Radio de Giro en Y (nm)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(paste(GRAPHDIR, "/gyr_y.svg", sep = ""), gyry_plot, width = 8, height = 6, dpi = 300)
# print(gyry_plot)

# Plot Z gyration radius
gyrz_plot <- ggplot(gyrbase, aes(x = Time)) +
  geom_line(aes(y = `Rg/sZ/N`, color = "Base"), linewidth = 0.5) +
  geom_line(data = gyrscaled, aes(y = `Rg/sZ/N`, color = "Rescaled"), linewidth = 0.5) +
  scale_color_manual(values = c("Base" = "blue", "Rescaled" = "#FF0000C0"),
                     labels = c("Base", "Reescalado"),
                     name = "Series") +
  scale_y_continuous(limits = c(1, 6), breaks = seq(1, 6, 1)) +
  labs(title = "Radio de Giro (eje Z) durante Producción",
       x = "Tiempo (ns)",
       y = "Radio de Giro en Z (nm)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(paste(GRAPHDIR, "/gyr_z.svg", sep = ""), gyrz_plot, width = 8, height = 6, dpi = 300)
#print(gyrz_plot)

### RMSD PLOT ###
rmsd_base <- readXVG(paste(BASEDIR, "/rmsd.xvg", sep = ""))
colnames(rmsd_base) <- c("Time", "RMSD Base")
rmsd_base$Time <- as.numeric(rmsd_base$Time)
rmsd_base$`RMSD Base` <- as.numeric(rmsd_base$`RMSD Base`)

rmsd_rescaled <- readXVG(paste(RESCALDIR, "/rmsd.xvg", sep = ""))
colnames(rmsd_rescaled) <- c("Time", "RMSD Rescaled")
rmsd_rescaled$Time <- as.numeric(rmsd_rescaled$Time)
rmsd_rescaled$`RMSD Rescaled` <- as.numeric(rmsd_rescaled$`RMSD Rescaled`)

rmsd <- merge(rmsd_base, rmsd_rescaled, by = "Time")
rm(rmsd_base, rmsd_rescaled) # Remove the data frames to free memory

# Transform time from ps to ns
rmsd$Time <- rmsd$Time / 1000

# Calculate mean and sd for the base and rescaled RMSD
mean_base <- signif(mean(rmsd$`RMSD Base`), 3)
mean_rescaled <- signif(mean(rmsd$`RMSD Rescaled`), 3)
sd_base <- signif(sd(rmsd$`RMSD Base`), 3)
sd_rescaled <- signif(sd(rmsd$`RMSD Rescaled`), 3)

# Print the mean and sd for the base and rescaled RMSD
cat("RMSD (nm):\n")
cat("Base: ", mean_base, " ± ", sd_base, "\n")
cat("Reescalado: ", mean_rescaled, " ± ", sd_rescaled, "\n")
cat("\n")

# Plot RMSD
rmsd_plot <- ggplot(rmsd, aes(x = Time)) +
  geom_line(aes(y = `RMSD Base`, color = "Base"), linewidth = 0.5) +
  geom_line(aes(y = `RMSD Rescaled`, color = "Rescaled"), linewidth = 0.5) +
  scale_color_manual(values = c("Base" = "blue", "Rescaled" = "#FF0000C0"),
                     labels = c("Base", "Reescalado"),
                     name = "Series") +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1)) +
  labs(title = "RMSD durante Producción",
       x = "Tiempo (ns)",
       y = "RMSD (nm)") +
  theme_minimal() +
  theme(legend.position = "right")

rmsd_plot_100ns <- rmsd_plot + scale_x_continuous(limits = c(0, 100)) +
  labs(title = "RMSD durante Producción (100 ns)")

ggsave(paste(GRAPHDIR, "/rmsd.svg", sep = ""), rmsd_plot, width = 8, height = 6, dpi = 300)
ggsave(paste(GRAPHDIR, "/rmsd_100ns.svg", sep = ""), rmsd_plot_100ns, width = 8, height = 6, dpi = 300)
# print(rmsd_plot)
# print(rmsd_plot_100ns)

### RMSF PLOT ###
rmsf_base <- readXVG(paste(BASEDIR, "/rmsf.xvg", sep = ""))
colnames(rmsf_base) <- c("Residuo", "RMSF Base")
rmsf_base$Residuo <- as.numeric(rmsf_base$Residuo)
rmsf_base$`RMSF Base` <- as.numeric(rmsf_base$`RMSF Base`)

rmsf_rescaled <- readXVG(paste(RESCALDIR, "/rmsf.xvg", sep = ""))
colnames(rmsf_rescaled) <- c("Residuo", "RMSF Rescaled")
rmsf_rescaled$Residuo <- as.numeric(rmsf_rescaled$Residuo)
rmsf_rescaled$`RMSF Rescaled` <- as.numeric(rmsf_rescaled$`RMSF Rescaled`)

# Fix the residue numbering, so it's the same as the row number
rmsf_base$Residuo <- 1:nrow(rmsf_base)
rmsf_rescaled$Residuo <- 1:nrow(rmsf_rescaled)

# Merge the data frames and remove the original ones
rmsf <- merge(rmsf_base, rmsf_rescaled, by = "Residuo")
rm(rmsf_base, rmsf_rescaled)

# Plot RMSF
rmsf_plot <- ggplot(rmsf, aes(x = Residuo)) +
  geom_line(aes(y = `RMSF Base`, color = "Base"), linewidth = 0.5) +
  geom_line(aes(y = `RMSF Rescaled`, color = "Rescaled"), linewidth = 0.5) +
  scale_color_manual(values = c("Base" = "blue", "Rescaled" = "#FF0000C0"),
                     labels = c("Base", "Reescalado"),
                     name = "Series") +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 0.5)) +
  labs(title = "RMSF durante Producción",
       x = "Residuo",
       y = "RMSF (nm)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(paste(GRAPHDIR, "/rmsf.svg", sep = ""), rmsf_plot, width = 8, height = 6, dpi = 300)
#print(rmsf_plot)