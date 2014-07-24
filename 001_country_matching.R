cou <- c("Algeria", "Italy")
affil <- c("1 place in, Algeria.", "2 place in Italy.", "3 place in, Ireland.", "4 Place in, Algeria.")

#returns true/false matrix if cou in affil
sapply(cou, grepl, affil)

#returns true/false lists for each element in affil
lapply(cou, grepl, affil)

#definitely loops through the list of affs, and gives true/false for every country
country.test <- lapply(query.affiliation, function(x) sapply(country.list, grepl, x))

#returns only true values, but not the value
#sub(pattern, replacement, text)
country.test <- sapply(query.affiliation, function(x) sub(pattern = country.list, replacement = any(sapply(country.list, grepl, x))))

#places counrty into a list
test1 <- if(grepl(cou[1], affil[1]) == TRUE) sub(affil[1], cou, cou[1])

