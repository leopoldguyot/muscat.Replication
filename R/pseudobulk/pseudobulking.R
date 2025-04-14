aggregate_assay <- function(data, method) {
    switch(method,
           "None" = data,
           "Mean" = aggregateData(x = data, fun = mean),
           "Sum" = aggregateData(x = data, fun = sum))
}
