query.aff.search <- data.frame()
df <- for(i in query.affiliation){
        if(grepl("Italy", i) == 1) {
                query.aff.search <- rbind(query.aff.search, i)
        }
        #else print("K")#query.aff.search <- rbind(query.aff.search, "Not Italy")
        
}