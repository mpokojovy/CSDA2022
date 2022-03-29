setwd("")

ncores = 200L

for (i in 1:ncores) {
  file_name = paste("chunk.", i, ".R", sep = "")
  file = file.create(file_name)
  write(paste("data.chunk = ", i, "\n", sep = ""), file = file_name, append = FALSE)
  write("set.seed(data.chunk)", file = file_name, append = TRUE)
  write(paste("setwd(\"", getwd(), "\")", sep = ""), file = file_name, append = TRUE)
  write("source(\"main.R\")", file = file_name, append = TRUE)
}